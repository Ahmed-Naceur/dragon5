*DECK T16RRE
      SUBROUTINE T16RRE(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                  NMXSR ,IMIREG,VELMTR,B2CRI ,BRNIRR,FLXINT,
     >                  OVERV,RECXSV,RECXSM,RECTMP,RECSCA)
*
*----
*
*Purpose:
*  Read tape16 REGION cross sections at a specific burnup.
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
* IMIREG  mixture update identifier where
*          =0 do not update;
*          =-1 update using CELLAV information;
*          > 0 update using specified region number.
* VELMTR  velocity for main transport.
* B2CRI   critical bucklings.
* FLXINT  volume integrated fluxes.
* OVERV   1/V cross sections.
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
      INTEGER          IFT16,IPRINT,NGCCPO,NGMTR,NVXSR,NMXSR,IMIREG
      INTEGER          IFGMTR(NGCCPO)
      REAL             VELMTR(NGMTR),B2CRI(3),BRNIRR(3),
     >                 FLXINT(NGCCPO),OVERV(NGCCPO),
     >                 RECXSV(NGCCPO,NVXSR+NMXSR),
     >                 RECXSM(NGCCPO,NGCCPO,NMXSR),
     >                 RECTMP(NGMTR,4),RECSCA(NGMTR,NGMTR)    
*----
*  T16 PARAMETERS
*----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=2)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NID,NJD
*----
*  LOCAL VARIABLES
*  WSMEV FACTOR TO TRANSFORM MEV IN JOULES (WS)
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             WSMEV
      PARAMETER       (IOUT=6,NAMSBR='T16RRE',WSMEV=1.602189E-13)
      INTEGER          IREG,IGR,IGC,IGD,IGF,JGR,JGC,JGD,JGF,
     >                 NREGON
      REAL             VOLUME,BRNTMP(3),RTIME
*----
*  INITIALIZE CROSS SECTION VECTORS
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDRSET(RECXSV,NGCCPO*(NVXSR+NMXSR),0.0)
      CALL XDRSET(RECXSM,NGCCPO*NGCCPO*NMXSR,0.0)
*----
*  LOCATE NEXT REGION DIMENSIONS RECORD
*  AND READ NREGON
*----
      IOPT=0
      TKEY1(1)='REGION    '
      TKEY2(1)='DIMENSIONS'
      NKEY=1
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. 2 ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,NREGON
      TKEY1(2)='CELLAV    '
      TKEY2(2)='NGROUPS   '
      NKEY=2
      DO IREG=1,NREGON
*----
*  REGIONAL FLUX
*----
        TKEY1(1)='REGION    '
        TKEY2(1)='FLUX      '
        CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >              NBE   )
        IF( NBE .NE. 3+NGMTR ) CALL XABORT(NAMSBR//
     >  ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
        IF(IMIREG .EQ. IREG) THEN
          READ(IFT16) RKEY1,RKEY2,NBE,NID,NJD,VOLUME,
     >               (RECTMP(IGR,1),IGR=1,NGMTR)
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6100) TKEY2(1)
            WRITE(IOUT,6110) (RECTMP(IGR,1),IGR=1,NGMTR)
          ENDIF
*----
*  TREAT ALL CONDENSED GROUPS
*----
          TKEY1(1)='REGION    '
          TKEY2(1)='SIGMAS    '
          IGF=0
          DO IGC=1,NGCCPO
            IGD=IGF+1
            IGF=IFGMTR(IGC)
*----
*  FLUX AND 1/V CROSS SECTION CONDENSATION
*----
            DO IGR=IGD,IGF
              FLXINT(IGC)=FLXINT(IGC)+RECTMP(IGR,1)
              OVERV(IGC)=OVERV(IGC)+RECTMP(IGR,1)/VELMTR(IGR)
            ENDDO
            IF(FLXINT(IGC) .NE. 0.0) THEN
              OVERV(IGC)=OVERV(IGC)/FLXINT(IGC)
              DO IGR=IGD,IGF
                RECTMP(IGR,1)=RECTMP(IGR,1)/FLXINT(IGC)
              ENDDO
              FLXINT(IGC)=FLXINT(IGC)*VOLUME
            ENDIF
*----
*  LOOP OBER MTR GROUP ASSOCIATED WITH CPO GROUPS
*----
            DO IGR=IGD,IGF
*----
*  READ CROSS SECTIONS
*----
              CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                    NBE   )
              IF( NBE .NE. 4+NGMTR ) CALL XABORT(NAMSBR//
     >        ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
              READ(IFT16) RKEY1,RKEY2,NBE,
     >                    RECTMP(IGR,4),RECTMP(IGR,3),RECTMP(IGR,2),
     >                    (RECSCA(IGR,JGR),JGR=1,NGMTR)
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6101) TKEY2(1),IGR
                WRITE(IOUT,6110)
     >            RECTMP(IGR,4),RECTMP(IGR,3),RECTMP(IGR,2),
     >           (RECSCA(IGR,JGR),JGR=1,NGMTR)
              ENDIF
*----
*  ABSORPTION, NU-FISSION AND TRANSPORT SECTION CONDENSATION
*----
              RECXSV(IGC, 2)=RECXSV(IGC, 2)
     >                      +RECTMP(IGR,2)*RECTMP(IGR,1)
              RECXSV(IGC, 3)=RECXSV(IGC, 3)
     >                      +RECTMP(IGR,3)*RECTMP(IGR,1)
              RECXSV(IGC,15)=RECXSV(IGC,15)
     >                      +RECTMP(IGR,4)*RECTMP(IGR,1)
*----
*  SCATTERING SECTION CONDENSATION
*----
              JGF=0
              DO JGC=1,NGCCPO
                JGD=JGF+1
                JGF=IFGMTR(JGC)
                DO JGR=JGD,JGF
                  RECXSM(JGC,IGC,1)=RECXSM(JGC,IGC,1)
     >                             +RECSCA(IGR,JGR)*RECTMP(IGR,1)
                  RECXSV(IGC,21)=RECXSV(IGC,21)
     >                          +RECSCA(IGR,JGR)*RECTMP(IGR,1)
                ENDDO
              ENDDO
            ENDDO
*----
*  TOTAL AND TRANSPORT CORRECTION
*----
            RECXSV(IGC,1)=RECXSV(IGC,15)+RECXSV(IGC,21)
            RECXSV(IGC,2)=RECXSV(IGC,1)-RECXSV(IGC,2)
          ENDDO
          IF( NBE .EQ. 2*NGMTR ) THEN
            IF(IPRINT .GE. 100) THEN
              RECTMP(IGR,3)=RECTMP(IGR,2)
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
     >             *(B2CRI(1)*RECTMP(IGR,2)+B2CRI(2)*RECTMP(IGR,3))
                RECXSV(IGC,18)=RECXSV(IGC,18)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,2)
                RECXSV(IGC,19)=RECXSV(IGC,19)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,2)
                RECXSV(IGC,20)=RECXSV(IGC,20)
     >                        +RECTMP(IGR,1)*RECTMP(IGR,3)
              ENDDO
              IF(RECXSV(IGC,17) .EQ. 0.0 .OR.
     >           RECXSV(IGC,18) .EQ. 0.0 .OR.
     >           RECXSV(IGC,19) .EQ. 0.0 .OR.
     >           RECXSV(IGC,19) .EQ. 0.0 ) THEN
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
          ELSE
            DO IGC=1,NGCCPO
              RECXSV(IGC,17)=1.0/(3.0*(RECXSV(IGC,1)-RECXSV(IGC,2)))
              RECXSV(IGC,18)=RECXSV(IGC,17)
              RECXSV(IGC,19)=RECXSV(IGC,17)
              RECXSV(IGC,20)=RECXSV(IGC,17)
            ENDDO
          ENDIF
          GO TO 105
        ELSE
          READ(IFT16) RKEY1,RKEY2,NBE
        ENDIF
      ENDDO
 105  CONTINUE
*----
*  READ FISSION SPECTRUM
*----
      TKEY1(1)='CELLAV    '
      TKEY2(1)='FISSPECT  '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NGMTR ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(RECTMP(IGR,4),IGR=1,NGMTR)
*----
*  CONDENSE FISSION SPECTRUM OVER CPO GROUPS
*----
      IGF=0
      DO IGC=1,NGCCPO
        IGD=IGF+1
        IGF=IFGMTR(IGC)
        DO IGR=IGD,IGF
          RECXSV(IGC, 5)=RECXSV(IGC,5)+RECTMP(IGR,4)
        ENDDO
      ENDDO
*----
*  BURNUP INFORMATION
*----
      TKEY1(2)='MTR       '
      TKEY2(2)='FEWGROUPS '
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
 6101 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP ',A10,
     >       6X,'GROUP  =',I10)
 6110 FORMAT(1P,10(2X,E10.3))
      END
