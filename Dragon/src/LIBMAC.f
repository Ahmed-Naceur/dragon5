*DECK LIBMAC
      SUBROUTINE LIBMAC(IPLIB ,IPLIBX,IPRINT,MAXISO,NBISO ,NBISOX,
     >                  IBSTEP,NBMIX ,NBMIXX,NGRO  ,TMPDAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read the information related to microscopic cross section library.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input/output
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPLIBX  pointer to the lattice old library density.
* IPRINT  print flag.
* MAXISO  maximum number of isotopes permitted.
* NBISO   number of isotopes present on IPLIB.
* NBISOX  number of isotopes present on IPLIBX.
* IBSTEP  burnup step on IPLIBX if 'BURN' option activated.
* NBMIX   number of mixtures defined on IPLIB.
* NBMIXX  number of mixtures defined on IPLIBX.
* NGRO    number of energy groups.
* TMPDAY  time/burnup/irradiation stamp in days.
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
      INTEGER      IOUT,NTC
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NTC=3,NAMSBR='LIBMAC')
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPLIB,IPLIBX
      INTEGER      IPRINT,MAXISO,NBISO,NBISOX,IBSTEP,NBMIX,NBMIXX,NGRO
      REAL         TMPDAY(3)
*----
*  INPUT
*----
      INTEGER      NSTATE
      PARAMETER   (NSTATE=40)
      INTEGER      ITYPLU,INTLIR,ILONG,ITYLCM,NCOMB,ISOT,IBM,J,
     >             ISTATE(NSTATE)
      CHARACTER    TEXT4*4,CARLIR*12
      REAL         REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL PARAMETERS
*----
      INTEGER      KCHAR(NTC),ISO,JSO,IMIX,NISOM,ITSTMP,NNMIX,MODISO,
     >             NMIXUP,NIUPD,IMIXX,ILCMLN,ILCMTY,ITC,ITEXT4
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,MIXIX,LOCUPD,LISM,IEVOL,
     > IEVOLX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,ISONMX
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,DENMIX,DENIX,DENMOD
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKM,MASKG
*----
*  SCRATCH STORAGE ALLOCATION
*   ISONAM  old name of isotopes.
*   ISONMX  new name of isotopes.
*   MIX     mix number of each isotope on IPLIB (can be zero).
*   DEN     density of each isotope on IPLIB.
*   DENMIX  density of mixture on IPLIB (can be -1.0).
*   MIXIX   mix number of each isotope on IPLIBX (can be zero).
*   DENIX   density of each isotope on IPLIBX.
*   LOCUPD  location of IPLIB mixture in IPLIBX.
*   LISM    location in IPLIB of isotope associated with a mixture.
*   DENMOD  modified density of each isotope on IPLIB.
*   MASKM   mixture update mask.
*   MASKG   group update mask.
*   IEVOL   flag making an isotope non-depleting:
*           =1 to force an isotope to be non-depleting;
*           =2 to force an isotope to be depleting;
*           =3 to force an isotope to be at saturation
*----
      ALLOCATE(MIX(MAXISO),MIXIX(NBISOX),LOCUPD(NBMIX),LISM(MAXISO),
     > IEVOL(MAXISO),IEVOLX(NBISOX))
      ALLOCATE(ISONAM(NTC,MAXISO),ISONMX(NTC,NBISOX))
      ALLOCATE(DEN(MAXISO),DENMIX(NBMIX),DENIX(NBISOX),DENMOD(MAXISO))
      ALLOCATE(MASKM(NBMIX),MASKG(NGRO))
*----
*  INITIALIZE
*----
      TEXT4='    '
      READ(TEXT4,'(A4)') ITEXT4
      CALL XDISET(ISONAM,NTC*MAXISO,ITEXT4)
      CALL XDISET(ISONMX,NTC*NBISOX,ITEXT4)
      CALL XDISET(MIX,MAXISO,0)
      CALL XDRSET(DEN,MAXISO,0.0)
      CALL XDRSET(DENMIX,NBMIX,-1.0)
*----
*  READ ORIGINAL ISOTOPE AND MIXTURE INFORMATION FROM IPLIB
*----
      IF(NBISO.GT.MAXISO) CALL XABORT('LIBMAC: MAXISO OVERFLOW.')
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DEN)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      CALL LCMLEN(IPLIB,'ISOTOPESTODO',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         CALL LCMGET(IPLIB,'ISOTOPESTODO',IEVOL)
      ELSE
         CALL XDISET(IEVOL,MAXISO,0)
      ENDIF
      DO 10 JSO=1,NBISO
        DENMOD(JSO)=DEN(JSO)
 10   CONTINUE
      CALL XDRSET(DENMIX,NBMIX,1.0)
      CALL XDISET(LOCUPD,NBMIX,0)
*----
*  WRITE ORIGINAL MATERIAL COMPOSITION IF REQUIRED
*----
      IF(IPRINT.GT.0) THEN
        WRITE(IOUT,6000) NGRO,NBISO,NBMIX
        DO 600 IMIX=1,NBMIX
          NISOM=0
          DO 601 ISO=1,NBISO
            IF(MIX(ISO).EQ.IMIX) THEN
               NISOM=NISOM+1
               LISM(NISOM)=ISO
            ENDIF
 601      CONTINUE
          IF(NISOM.GT.0) THEN
            WRITE(IOUT,6010) IMIX
            WRITE(IOUT,6011) ((ISONAM(ITC,LISM(ISO)),ITC=1,NTC-1),
     >        DEN(LISM(ISO)),IEVOL(LISM(ISO)),ISO=1,NISOM)
          ENDIF
 600    CONTINUE
      ENDIF
*----
*  READ ISOTOPE AND MIXTURE INFORMATION FROM IPLIBX
*----
      ITSTMP=2
      IF((IBSTEP.EQ.0).AND.(NBISOX.GT.0)) THEN
*       READ FROM A MICROLIB.
        CALL LCMGET(IPLIBX,'ISOTOPESDENS',DENIX)
        CALL LCMGET(IPLIBX,'ISOTOPESMIX',MIXIX)
        CALL LCMGET(IPLIBX,'ISOTOPESUSED',ISONMX)
        CALL LCMLEN(IPLIBX,'ISOTOPESTODO',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
           CALL LCMGET(IPLIBX,'ISOTOPESTODO',IEVOLX)
        ELSE
           CALL XDISET(IEVOLX,NBISOX,0)
        ENDIF
      ELSE IF((IBSTEP.GT.0).AND.(NBISOX.GT.0)) THEN
*       READ FROM A BURNUP OBJECT.
        WRITE(CARLIR,'(8HDEPL-DAT,I4.4)') IBSTEP
        CALL LCMGET(IPLIBX,'ISOTOPESMIX',MIXIX)
        CALL LCMGET(IPLIBX,'ISOTOPESUSED',ISONMX)
        CALL LCMLEN(IPLIBX,'ISOTOPESTODO',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
           CALL LCMGET(IPLIBX,'ISOTOPESTODO',IEVOLX)
        ELSE
           CALL XDISET(IEVOLX,NBISOX,0)
        ENDIF
        CALL LCMSIX(IPLIBX,CARLIR,1)
        CALL LCMGET(IPLIBX,'ISOTOPESDENS',DENIX)
        CALL LCMLEN(IPLIBX,'BURNUP-IRRAD',ILONG,ITYLCM)
        IF(ILONG.EQ.2) THEN
           CALL LCMGET(IPLIBX,'BURNUP-IRRAD',TMPDAY(2))
        ENDIF
        CALL LCMSIX(IPLIBX,' ',2)
      ENDIF
      IF(IPRINT.GT.0) THEN
        WRITE(IOUT,6001) NGRO,NBISOX,NBMIX
        DO 620 IMIX=1,NBMIX
          NISOM=0
          DO 621 ISO=1,NBISO
            IF(MIX(ISO).EQ.IMIX) THEN
               NISOM=NISOM+1
               LISM(NISOM)=ISO
            ENDIF
 621      CONTINUE
          WRITE(IOUT,6010) IMIX
          WRITE(IOUT,6011) ((ISONAM(ITC,LISM(ISO)),ITC=1,NTC-1),
     >      DENMOD(LISM(ISO)),IEVOLX(LISM(ISO)),ISO=1,NISOM)
 620    CONTINUE
      ENDIF
*----
*  READ UPDATE INFORMATION FROM INPUT
*  FORMAT PERMITTED ARE
*  [MIX IMLIB [ IMLIBX ] [ DENMOD ] [ NAMISO CONCM(I) ] ] ;
*  DEFAULT:
*    MIX ABSENT, IPLIBX >0     -> ALL ISOTOPES AND ALL MIXTURES
*    MIX ABSENT, IPLIBX =0     -> NO UPDATE - PRINT ONLY
*    IMLIBX ABSENT, IPLIBX >0  -> CONCF(ISO,IMX)=CONC(ISO,IMLIBX)
*    IMLIBX ABSENT, IPLIBX =0  -> CONCF(ISO,IMX)=CONC(ISO,-IMLIB)
*    IMLIBX > 0                -> CORRECTION FROM IPLIBX
*    IMLIBX = -IMLIB           -> CORRECTION FROM IPLIB
*    DENMOD PRESENT            -> CONCF(I,IMX)=CONCF(I,IMX)*DENMOD
*    NAMISO ABSENT, IPLIBX >0  -> ALL ISOTOPE FOR MIXTURE
*    NAMISO ABSENT, IPLIBX =0  -> ALL ISOTOPE FOR MIXTURE
*    NAMISO PRESENT, CONCM >=0 -> ISOTOPE SPECIFIED
*                                 CONCF(I,IMX)=DENMOD*CONCM(I,IMX)
*    NAMISO PRESENT, CONCM <0  -> ISOTOPE SPECIFIED
*                                 CONCF(I,IMX)=CONC(ISO,-IMLIB)*DENMOD
*----
      NNMIX=0
 100  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER DATA EXPECTED.')
      IF(CARLIR.EQ.';') THEN
        GO TO 105
      ELSE IF(CARLIR(1:3).EQ.'MIX') THEN
        CALL REDGET(ITYPLU,NNMIX,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': MIXTURE TO UPDATE MUST BE GIVEN.')
        IF(NNMIX.GT.NBMIX) THEN
          CALL XABORT(NAMSBR//': CANNOT UPDATE THIS MIXTURE.')
        ELSE IF(NNMIX.LE.0) THEN
          CALL XABORT(NAMSBR//': MIX NUMBER.LE.0.')
        ENDIF
        IF(IBSTEP.EQ.0) THEN
          LOCUPD(NNMIX)=NNMIX
        ELSE
          LOCUPD(NNMIX)=-NNMIX
        ENDIF
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          IF(INTLIR.LE.0 .OR.
     >       INTLIR.GT.NBMIXX) CALL XABORT(NAMSBR//
     >       ': CANNOT UPDATE THIS MIXTURE.')
          LOCUPD(NNMIX)=INTLIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
        IF(ITYPLU.EQ.2) THEN
          DENMIX(NNMIX)=REALIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
        GO TO 101
      ENDIF
      IF(CARLIR.EQ.'NOEV') THEN
        DO 30 ISO=1,NBISO
        IF(NNMIX.EQ.MIX(ISO)) IEVOL(ISO)=1
 30     CONTINUE
      ELSE IF(CARLIR.EQ.'EVOL') THEN
        DO 35 ISO=1,NBISO
        IF(NNMIX.EQ.MIX(ISO)) IEVOL(ISO)=2
 35     CONTINUE
      ELSE
        READ(CARLIR,'(2A4)') (KCHAR(ITC),ITC=1,NTC-1)
        MODISO=0
        IF(LOCUPD(NNMIX).LT.0) THEN
          DO 40 ISO=1,NBISO
            IF(KCHAR(1).EQ.ISONAM(1,ISO) .AND.
     >         KCHAR(2).EQ.ISONAM(2,ISO) .AND.
     >         NNMIX.EQ.MIX(ISO)) THEN
               MODISO=ISO
               GO TO 45
            ENDIF
 40       CONTINUE
          WRITE(IOUT,'(10H MIXTURE :,1X,I10,10H ISOTOPE :,1X,2A4)')
     >      NNMIX,(KCHAR(ITC),ITC=1,NTC-1)
          CALL XABORT(NAMSBR//
     >    ': CANNOT UPDATE THIS ISOTOPE IN CURRENT MIXTURE.')
 45       CONTINUE
          CALL REDGET(ITYPLU,INTLIR,DENMOD(MODISO),CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': NEW ISOTOPIC DENSITY EXPECTED.')
        ELSE IF(LOCUPD(NNMIX).GT.0) THEN
          DO 50 ISO=1,NBISOX
            IF(KCHAR(1).EQ.ISONMX(1,ISO) .AND.
     >         KCHAR(2).EQ.ISONMX(2,ISO) .AND.
     >         LOCUPD(NNMIX).EQ.MIXIX(ISO)) THEN
             MODISO=ISO
             GO TO 55
          ENDIF
 50       CONTINUE
          WRITE(IOUT,'(10H MIXTURE :,1X,I10,10H ISOTOPE :,1X,2A4)')
     >      NNMIX,(KCHAR(ITC),ITC=1,NTC-1)
          CALL XABORT(NAMSBR//
     >    ': CANNOT UPDATE THIS ISOTOPE IN CURRENT MIXTURE.')
 55       CONTINUE
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': NEW ISOTOPIC DENSITY EXPECTED.')
          IF(REALIR.LT.0.0) THEN
            DENIX(MODISO)=-DENIX(MODISO)
          ELSE
            DENIX(MODISO)=REALIR
          ENDIF
        ENDIF
      ENDIF
      GO TO 100
 105  CONTINUE
      IF((NNMIX.EQ.0).AND.(NBISOX.GT.0)) THEN
        IF(IBSTEP.EQ.0) THEN
          CALL XDISET(LOCUPD,NBMIX,IMIX)
        ELSE IF(IBSTEP.GT.0) THEN
          CALL XDISET(LOCUPD,NBMIX,-IMIX)
        ENDIF
      ENDIF
*----
*  TRANSFER DENSITY FROM DENMOD OR DENIX TO DEN IF REQUIRED
*----
      NMIXUP=0
      DO 70 IMIX=1,NBMIX
        MASKM(IMIX)=.FALSE.
        IF(LOCUPD(IMIX).GT.0) THEN
          NIUPD=0
          IMIXX=LOCUPD(IMIX)
          DO 71 ISO=1,NBISOX
            IF(MIXIX(ISO).EQ.IMIXX) THEN
              DO 72 JSO=1,NBISO
                IF(ISONAM(1,JSO).EQ.ISONMX(1,ISO) .AND.
     >             ISONAM(2,JSO).EQ.ISONMX(2,ISO) .AND.
     >             MIX(JSO)     .EQ.IMIX) THEN
                  IF(DENMIX(IMIX)*DENIX(ISO) .NE. DEN(JSO)) THEN
                    DEN(JSO)=DENMIX(IMIX)*DENIX(ISO)
                    NIUPD=NIUPD+1
                  ENDIF
                ENDIF
 72           CONTINUE
            ENDIF
 71       CONTINUE
          IF(NIUPD .NE. 0) THEN
            MASKM(IMIX)=.TRUE.
            NMIXUP=NMIXUP+1
          ENDIF
        ELSE IF(LOCUPD(IMIX).LT.0) THEN
          NIUPD=0
          DO 73 ISO=1,NBISOX
            IF(MIXIX(ISO).EQ.IMIX) THEN
              DO 74 JSO=1,NBISO
                IF(ISONAM(1,JSO).EQ.ISONMX(1,ISO) .AND.
     >             ISONAM(2,JSO).EQ.ISONMX(2,ISO) .AND.
     >             MIX(JSO)     .EQ.IMIX) THEN
                  IF(DENMIX(IMIX)*DENIX(ISO) .NE. DEN(JSO)) THEN
                    IF(DENMIX(IMIX)*DENIX(ISO).GE.0.0) THEN
                      DEN(JSO)=DENMIX(IMIX)*DENIX(ISO)
                      NIUPD=NIUPD+1
                    ENDIF
                  ENDIF
                ENDIF
 74           CONTINUE
            ENDIF
 73       CONTINUE
          IF(NIUPD.NE.0) THEN
            MASKM(IMIX)=.TRUE.
            NMIXUP=NMIXUP+1
          ENDIF
        ENDIF
 70   CONTINUE
*----
*  UPDATE ALL MATERIAL IF MACROLIB DIRECTORY ABSENT
*----
      CALL LCMLEN(IPLIB,'MACROLIB',ILCMLN,ILCMTY)
      IF(ILCMLN.EQ.0) THEN
        NMIXUP=NBMIX
        CALL XDLSET(MASKM,NBMIX,.TRUE.)
        IF(IPRINT.GT.0) WRITE(IOUT,6004)
      ENDIF
*----
*  RECOMPUTE THE NUMBER OF DEPLETING MIXTURES
*----
      IF(NMIXUP.GT.0) THEN
        NCOMB=0
        DO 90 ISOT=1,NBISO
        IBM=MIX(ISOT)
        IF(IBM.EQ.0) GO TO 90
        IF(IEVOL(ISOT).NE.1) THEN
           DO 80 J=1,NCOMB
           IF(IBM.EQ.LOCUPD(J)) GO TO 90
   80      CONTINUE
           NCOMB=NCOMB+1
           LOCUPD(NCOMB)=IBM
           GO TO 90
        ENDIF
   90   CONTINUE
        IF(IPRINT.GT.0) WRITE(IOUT,6020) NCOMB
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        IF(ISTATE(12).NE.NCOMB) THEN
           ISTATE(12)=NCOMB
           CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        ENDIF
*----
*  WRITE UPDATED MATERIAL COMPOSITION IF REQUIRED
*----
        IF(IPRINT.GT.0) THEN
          WRITE(IOUT,6002)
          DO 630 IMIX=1,NBMIX
            NISOM=0
            IF(MASKM(IMIX)) THEN
              DO 631 ISO=1,NBISO
                IF(MIX(ISO).EQ.IMIX) THEN
                   NISOM=NISOM+1
                   LISM(NISOM)=ISO
                ENDIF
 631          CONTINUE
              IF(NISOM.GT.0) THEN
                WRITE(IOUT,6010) IMIX
                WRITE(IOUT,6011) ((ISONAM(ITC,LISM(ISO)),ITC=1,NTC-1),
     >              DEN(LISM(ISO)),IEVOL(LISM(ISO)),ISO=1,NISOM)
              ENDIF
            ENDIF
 630      CONTINUE
        ENDIF
*----
*  SAVE ISOTOPE NEW DENSITY
*----
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DEN)
        CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO,1,IEVOL)
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS
*----
        CALL XDLSET(MASKG,NGRO,.TRUE.)
        CALL LIBMIX(IPLIB,NBMIX,NGRO,NBISO,ISONAM,MIX,DEN,MASKM,MASKG,
     >  ITSTMP,TMPDAY)
      ELSE
        IF(IPRINT.GT.0) WRITE(IOUT,6003)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MASKG,MASKM)
      DEALLOCATE(DENMOD,DENIX,DENMIX,DEN)
      DEALLOCATE(ISONMX,ISONAM)
      DEALLOCATE(IEVOLX,IEVOL,LISM,LOCUPD,MIXIX,MIX)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/' LIBMAC: MODIFIED LIBRARY PROPERTIES '/
     >        '         NUMBER OF GROUPS   = ',I10/
     >        '         NUMBER OF ISOTOPES = ',I10/
     >        '         NUMBER OF MIXTURES = ',I10/
     >        '         ORIGINAL NUMBER DENSITIES IN MIXTURES',
     >        ' FOLLOWS')
 6001 FORMAT(/' LIBMAC: OLD LIBRARY PROPERTIES (READ ONLY) '/
     >        '         NUMBER OF GROUPS   = ',I10/
     >        '         NUMBER OF ISOTOPES = ',I10/
     >        '         NUMBER OF MIXTURES = ',I10/
     >        '         ORIGINAL NUMBER DENSITIES IN MIXTURES',
     >        ' FOLLOWS')
 6002 FORMAT(/' LIBMAC: FINAL NUMBER DENSITIES MIXTURES FOLLOWS')
 6003 FORMAT(/' LIBMAC: NO UPDATED MIXTURES')
 6004 FORMAT(/' LIBMAC: MACROSCOPIC ABSENT -> ALL MIXTURES UPDATED')
 6010 FORMAT(/' ISOTOPIC DENSITIES FOR MIXTURE =',I4)
 6011 FORMAT(1P,5(4X,2A4,':',E12.4,' (',I1,')'))
 6020 FORMAT(/' LIBMAC: NUMBER OF DEPLETING MIXTURES =',I4)
      END
