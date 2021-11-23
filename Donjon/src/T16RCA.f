*DECK T16RCA
      SUBROUTINE T16RCA(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                  NMXSR ,B2CRI ,BRNIRR, NZONE,RECXSV,RECXSM,
     >                  RECTMP,RECSCA,ZONVOL)
*
*----
*
*Purpose:
*  Read tape16 CELLAV cross sections at a specific burnup.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         =0 for no print; >=  1 print processing option.
* NGCCPO  number of edit groups.
* NGMTR   number of main transport groups.
* IFGMTR  fewgroups for main transport.
* NVXSR   number of vector cross sections.
* NMXSR   number of matrix cross sections.
* B2CRI   critical bucklings.
* NZONE   number of zones.
* ZONVOL  zone volume.
*
*Parameters: output
* BRNIRR  burnup and irradiation.
* RECXSV  vector cross sections.
* RECXSM  matrix cross sections.
* RECTMP  dummy vector cross sections.
* RECSCA  dummy matrix cross sections.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,NGCCPO,NGMTR,NVXSR,NMXSR
      INTEGER          IFGMTR(NGCCPO)
      REAL             B2CRI(3),BRNIRR(3),
     >                 RECXSV(NGCCPO,NVXSR+NMXSR),
     >                 RECXSM(NGCCPO,NGCCPO,NMXSR),
     >                 RECTMP(NGMTR,4),RECSCA(NGMTR,NGMTR)
      REAL             ZONVOL(NZONE)
*----
*  T16 PARAMETERS
*----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=2)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NID,IR, IZ
      REAL             RID
*----
*  LOCAL VARIABLES
*  WSMEV FACTOR TO TRANSFORM MEV IN JOULES (WS)
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             WSMEV
      PARAMETER       (IOUT=6,NAMSBR='T16RCA',WSMEV=1.602189E-13)
      INTEGER          IGR,IGC,IGD,IGF,JGR,JGC,JGD,JGF
      REAL             FLXNOR,BRNTMP(3),RTIME
      INTEGER          NZONE
      REAL             CELLV
*----
*  INITIALIZE CROSS SECTION VECTORS
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDRSET(RECXSV,NGCCPO*(NVXSR+NMXSR),0.0)
      CALL XDRSET(RECXSM,NGCCPO*NGCCPO*NMXSR,0.0)
*----
*  LOCATE NEXT CELLAV RECORD
*----
      IOPT=0
      TKEY1(1)='CELLAV    '
      TKEY2(1)='MODERATOR '
      NKEY=1
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .LE. 0 ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE
      NKEY=2
*----
*  CELL AVERAGED ABSORPTION X-S
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='ABSORPTION'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,4),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,4),IGR=1,NGMTR)
      ENDIF
*----
*  CELL AVERAGED NU*FISSION
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='NU-FISSION'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,3),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,3),IGR=1,NGMTR)
      ENDIF
*----
*  CELL AVERAGED TRANSPORT
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='TOTAL-X   '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,2),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,2),IGR=1,NGMTR)
      ENDIF
      CELLV=0.0
      DO IZ=1, NZONE 
        CELLV=CELLV+ZONVOL(IZ)
      ENDDO 
      TKEY1(1)='CELLAV    '
      TKEY2(1)='FLUX   '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,1),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,1),IGR=1,NGMTR)
      ENDIF
      DO IGR=1, NGMTR 
        RECTMP(IGR,1)=RECTMP(IGR,1)/CELLV
      ENDDO 
*----
*  CONDENSE TRANSPORT, ABSORPTION AND NU-FISSION X-S
*  OVER CPO GROUPS
*----
      IGF=0
      DO IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        FLXNOR=0.0
        DO IGR=IGD,IGF
          FLXNOR=FLXNOR+RECTMP(IGR,1)
        ENDDO
        IF(FLXNOR .GT. 0.0) THEN
          FLXNOR=1.0/FLXNOR
          DO IGR=IGD,IGF
            RECTMP(IGR, 1)=RECTMP(IGR, 1)*FLXNOR
            RECXSV(IGC, 2)=RECXSV(IGC, 2)
     >                    +RECTMP(IGR,2)*RECTMP(IGR,1)
            RECXSV(IGC, 3)=RECXSV(IGC, 3)
     >                    +RECTMP(IGR,3)*RECTMP(IGR,1)
            RECXSV(IGC,15)=RECXSV(IGC,15)
     >                    +RECTMP(IGR,4)*RECTMP(IGR,1)
          ENDDO
        ELSE
          CALL XABORT(NAMSBR//
     >   ': FLUX IN ONE CPO GROUP IS 0.0')
        ENDIF
      ENDDO
*----
*  ISOTROPIC SCATTERING MATRIX FROM GROUP IGR TO JGR
*  IS STORED ON TAPE 16 AS
*  ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
*  RECXSM(IGTO,IGFROM,1) REPRESENT
*  SCATTERING CROSS SECTION
*  FROM GROUP "IGFROM" TO GROUP "IGTO"
*  FOR ANISOTROPY LEVEL 1
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='SCATTER   '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR*NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,
     >          ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
      ENDIF
*----
*  FISSION SPECTRUM
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='FISSPECT  '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,4),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,4),IGR=1,NGMTR)
      ENDIF
*----
*  CONDENSE ISOTROPIC SCATTERING MATRIX AND FISSION SPECTRUM
*  OVER CPO GROUPS
*  COMPUTE TOTAL ISOTROPIC SCATTERING
*  COMPUTE TOTAL AND TRANSPORT CORRECTION
*  TOTAL(1) = ABSORPTION (15) + SCATTERING (21)
*  TRANSPORT CORRECTION (2) = TOTAL(1) -TRANSPORT CORRECTED (2)
*----
      IGF=0
      DO IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        DO IGR=IGD,IGF
          RECXSV(IGC, 5)=RECXSV(IGC,5)+RECTMP(IGR,4)
          JGF=0
          DO JGC=1,NGCCPO
            JGD=JGF+1
            JGF=IFGMTR(JGC)
            DO JGR=JGD,JGF
              RECXSM(JGC,IGC,1)=RECXSM(JGC,IGC,1)
     >                         +RECSCA(IGR,JGR)*RECTMP(IGR,1)
              RECXSV(IGC,21)=RECXSV(IGC,21)
     >                      +RECSCA(IGR,JGR)*RECTMP(IGR,1)
            ENDDO
          ENDDO
        ENDDO
        RECXSV(IGC,1)=RECXSV(IGC,15)+RECXSV(IGC,21)
        RECXSV(IGC,2)=RECXSV(IGC,1)-RECXSV(IGC,2)
      ENDDO
*----
*  LINEARLY ANISOTROPIC SCATTERING FROM GROUP IGR TO JGR
*  IS STORED ON TAPE 16 AS
*  ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
*  RECXSM(IGTO,IGFROM,2) REPRESENT
*  SCATTERING CROSS SECTION
*  FROM GROUP "IGFROM" TO GROUP "IGTO"
*  FOR ANISOTROPY LEVEL 2
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='SCATERP1  '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .EQ. NGMTR*NGMTR ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,
     >            ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6100) TKEY2(1)
          WRITE(IOUT,6110) ((RECSCA(IGR,JGR),IGR=1,NGMTR),JGR=1,NGMTR)
        ENDIF
*----
*  CONDENSE LINEARLY ANISOTROPIC SCATTERING MATRIX
*  OVER CPO GROUPS
*  COMPUTE TOTAL LINEARLY ANISOTROPIC SCATTERING
*----
        IGF=0
        DO IGC=1,NGCCPO
          IGD=IGF+1
          IGF=IFGMTR(IGC)
          DO IGR=IGD,IGF
            JGF=0
            DO JGC=1,NGCCPO
              JGD=JGF+1
              JGF=IFGMTR(JGC)
              DO JGR=JGD,JGF
                RECXSM(JGC,IGC,2)=RECXSM(JGC,IGC,2)
     >                           +RECTMP(IGR,4)*RECSCA(IGR,JGR)
                RECXSV(IGC,22)=RECXSV(IGC,22)
     >                        +RECTMP(IGR,4)*RECSCA(IGR,JGR)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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
     >           (RECTMP(IGR,2),IGR=1,NGMTR),
     >           (RECTMP(IGR,3),IGR=1,NGMTR),
     >           (RID,IGR=1,NGMTR),
     >           (RID,IR=1,2)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,2),IGR=1,NGMTR)
        WRITE(IOUT,6110) (RECTMP(IGR,3),IGR=1,NGMTR)
      ENDIF
*----
*  CONDENSE DIFFUSION COEFFICIENTS
*  COMPUTE STRD=1/3*DIFF
*----
      IGF=0
      DO IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        DO IGR=IGD,IGF
          RECXSV(IGC,17)=RECXSV(IGC,17)+RECTMP(IGR,1)
     >       *(B2CRI(1)*RECTMP(IGR,2)+B2CRI(2)*RECTMP(IGR,3))
          RECXSV(IGC,18)=RECXSV(IGC,18)
     >                  +RECTMP(IGR,1)*RECTMP(IGR,2)
          RECXSV(IGC,19)=RECXSV(IGC,19)
     >                  +RECTMP(IGR,1)*RECTMP(IGR,2)
          RECXSV(IGC,20)=RECXSV(IGC,20)
     >                  +RECTMP(IGR,1)*RECTMP(IGR,3)
        ENDDO
*----
*  IF DIFFUSION COEFFICIENT VANISHES
*  ASSUME D=1/3*(TRANSPORT CORRECTED)
*  NO DIRECTIONAL EFFECT
*  THEN  USE STRD=1/3*DIFF
*----
        IF(RECXSV(IGC,17) .EQ. 0.0 .OR.
     >     RECXSV(IGC,18) .EQ. 0.0 .OR.
     >     RECXSV(IGC,19) .EQ. 0.0 .OR.
     >     RECXSV(IGC,19) .EQ. 0.0 ) THEN
          RECXSV(IGC,17)=RECXSV(IGC,1)-RECXSV(IGC,2)
          RECXSV(IGC,18)=0.0
          RECXSV(IGC,19)=0.0
          RECXSV(IGC,20)=0.0
        ELSE
          RECXSV(IGC,17)=1.0/(3.0*RECXSV(IGC,17))
          RECXSV(IGC,18)=1.0/(3.0*RECXSV(IGC,18))
          RECXSV(IGC,19)=1.0/(3.0*RECXSV(IGC,19))
          RECXSV(IGC,20)=1.0/(3.0*RECXSV(IGC,20))
        ENDIF
      ENDDO
*----
*  FISSION CROSS SECTION
*----
      TKEY1(2)='MTR       '
      TKEY2(2)='FEWGROUPS '
      TKEY1(1)='CELLAV    '
      TKEY2(1)='SIGMAF    '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,4),IGR=1,NGMTR)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6100) TKEY2(1)
        WRITE(IOUT,6110) (RECTMP(IGR,4),IGR=1,NGMTR)
      ENDIF
*----
*  CONDENSE FISSION CROSS SECTION
*  OVER CPO GROUPS
*----
      IGF=0
      DO IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        DO IGR=IGD,IGF
          RECXSV(IGC, 4)=RECXSV(IGC, 4)
     >                  +RECTMP(IGR,4)*RECTMP(IGR,1)
        ENDDO
      ENDDO
*----
*  BURNUP INFORMATION
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='AVG-ENERGY'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .EQ. 5 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,RTIME,
     >              BRNTMP(3),BRNTMP(1),BRNTMP(2)
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6010) RTIME,BRNTMP(3),BRNTMP(1),BRNTMP(2)
        ENDIF
        BRNIRR(1)=BRNTMP(1)
        BRNIRR(2)=BRNTMP(2)
        BRNIRR(3)=WSMEV*BRNTMP(3)
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001)
      ENDIF
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*'))
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'BURNUP IRRADIATION '/1P,
     >       6X,'TIME    (DAYS)     = ',E10.3/
     >       6X,'ENERGY  (MEV)      = ',E10.3/
     >       6X,'BURNUP  (MWD/T)    = ',E10.3/
     >       6X,'IRRADIATION (N/KB) = ',E10.3)
 6100 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP ',A10)
 6110 FORMAT(1P,10(2X,E10.3))
      END
