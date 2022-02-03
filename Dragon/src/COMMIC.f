*DECK COMMIC
      SUBROUTINE COMMIC(IMPX,IPCPO,IPEDIT,IPEDI2,LMACRO,ICAL,MAXCAL,
     1 NMIL,NISOTS,NG,NED,NW,FNORM,LISO,NISOP,NOMISP,NGFF,NALBP,IDF,
     2 ITRES)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover a microlib corresponding to a set of homogenized mixtures.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* IPCPO   pointer to the multicompo.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* IPEDI2  pointer to the edition object containing group form factor
*         information (L_EDIT signature).
* LMACRO  flag set to .TRUE. to recover cross sections from the
*         macrolib.
* ICAL    index of the elementary calculation.
* MAXCAL  maximum number of elementary calculations in the multicompo.
* NMIL    number of homogenized mixtures.
* NISOTS  number of isotopes in the microlib pointed by IPEDIT.
* NG      number of energy groups.
* NED     number of additional edits.
* NW      type of weighting for P1 cross section info (=0: P0 ; =1: P1).
* FNORM   flux normalization factor.
* LISO    =.true. if we want to register the region number of the
*         isotopes.
* NISOP   number of user-requested particularized isotopes. Equal to
*         zero if all EDI: isotopes are particularized.
* NOMISP  names of user-requested particularized isotopes.
* NGFF    number of form factors per energy group.
* NALBP   number of physical albedos per energy group.
* IDF     flag for ADF info (-1/0/1/2: candidate/absent/present).
*
*Parameters: output
* ITRES   creation index for the macroscopic residual (=0: not created;
*         =1: not a FP precursor; =2: is a FP precursor).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO,IPEDIT,IPEDI2
      INTEGER IMPX,ICAL,MAXCAL,NMIL,NISOTS,NG,NED,NW,NISOP,NGFF,NALBP,
     1 IDF,ITRES
      CHARACTER NOMISP(NISOP)*8
      REAL FNORM
      LOGICAL LMACRO,LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO,MPCPO,NPCPO,OPCPO,IPWORK,JPEDIT,
     1 KPEDIT
      CHARACTER TEXT4*4,TEXT8*8,TEXT12*12
      INTEGER IPAR(NSTATE),ISTATE(NSTATE)
      LOGICAL LRES
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX1,ITYP1,ITOD1,ITYP2,
     1 ITOD2,IPIFI,IPIFI2,ISW
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HUSE1,HNAM1,HUSE2,HNAM2,
     1 HVECT
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS1,TEMP1,VOL1,DENS2,TEMP2,
     1 VOL2,WORK,ENER,ZLAMB,DELT,VOLMIX,PYIELD,PYIEL2,PYRES
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ALBP,ADF,ADF2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ALBP2
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(HUSE1(3,NISOTS),HNAM1(3,NISOTS),IMIX1(NISOTS),
     1 ITYP1(NISOTS),ITOD1(NISOTS),HUSE2(3,NISOTS),HNAM2(3,NISOTS),
     2 ITYP2(NISOTS),ITOD2(NISOTS),ISW(NISOTS),HVECT(2,NED+1))
      ALLOCATE(DENS1(NISOTS),TEMP1(NISOTS),VOL1(NISOTS),DENS2(NISOTS),
     1 TEMP2(NISOTS),VOL2(NISOTS),WORK(NG),ENER(NG+1),DELT(NG),
     2 VOLMIX(NMIL))
*----
*  RECOVER THE TOC RECORDS
*----
      IF(.NOT.LMACRO) THEN
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
         CALL LCMGET(IPEDIT,'ISOTOPESUSED',HUSE1)
         CALL LCMGET(IPEDIT,'ISOTOPERNAME',HNAM1)
         CALL LCMGET(IPEDIT,'ISOTOPESDENS',DENS1)
         CALL LCMGET(IPEDIT,'ISOTOPESMIX',IMIX1)
         CALL LCMGET(IPEDIT,'ISOTOPESTYPE',ITYP1)
         CALL LCMGET(IPEDIT,'ISOTOPESTODO',ITOD1)
         CALL LCMGET(IPEDIT,'ISOTOPESVOL',VOL1)
         CALL LCMGET(IPEDIT,'ISOTOPESTEMP',TEMP1)
         CALL LCMGET(IPEDIT,'MIXTURESVOL',VOLMIX)
         CALL LCMGET(IPEDIT,'K-EFFECTIVE',EIGENK)
         CALL LCMGET(IPEDIT,'K-INFINITY',EIGINF)
         CALL LCMLEN(IPEDIT,'B2  B1HOM',ILONG,ITYLCM)
         IF(ILONG.EQ.1) THEN
            CALL LCMGET(IPEDIT,'B2  B1HOM',B2)
         ELSE
            B2=0.0
         ENDIF
         IF(NED.GT.0) CALL LCMGET(IPEDIT,'ADDXSNAME-P0',HVECT)
         CALL LCMGET(IPEDIT,'ENERGY',ENER)
         CALL LCMGET(IPEDIT,'DELTAU',DELT)
         NDEL=IPAR(19)
         IF(IDF.NE.0) IDF=IPAR(24)
      ENDIF
*----
*  LOOP OVER HOMOGENIZED MIXTURES
*----
      TEXT4='    '
      NSPH=0
      READ(TEXT4,'(A4)') ITEXT
      JPCPO=LCMLID(IPCPO,'MIXTURES',NMIL)
      DO 130 IMIL=1,NMIL
      KPCPO=LCMDIL(JPCPO,IMIL)
      LPCPO=LCMLID(KPCPO,'CALCULATIONS',MAXCAL)
      MPCPO=LCMDIL(LPCPO,ICAL)
      ISO3=0
      NPCPO=LCMLID(MPCPO,'ISOTOPESLIST',NISOTS)
      IF(LMACRO) THEN
*        RECOVER CROSS SECTIONS FROM THE MACROLIB.
         CALL LCMSIX(IPEDIT,'MACROLIB',1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
         NL=IPAR(3)
         NF=IPAR(4)
         NDEL=IPAR(7)
         IF(IDF.NE.0) IDF=IPAR(12)
         NSPH=IPAR(14)
         ALLOCATE(ZLAMB(NDEL))
         OPCPO=LCMDIL(NPCPO,1) ! set first isotope
         CALL COMACR(IPEDIT,IMPX,OPCPO,NG,NMIL,NED,NL,NF,NDEL,NW,IMIL,
     1   FNORM,NSPH,EIGENK,EIGINF,B2,VOLUME,ENER,DELT,HVECT,ZLAMB)
         CALL LCMSIX(IPEDIT,' ',2)
         NISO2=1
         DENS2(NISO2)=1.0
         ITYP2(NISO2)=1
         IF(NF.GT.0) ITYP2(NISO2)=2
         ITOD2(NISO2)=1
         VOL2(NISO2)=VOLUME
         VOLMIX(IMIL)=VOLUME
         TEMP2(NISO2)=0.0
         TEXT12='*MAC*RES'
         CALL LCMPTC(OPCPO,'ALIAS',12,1,TEXT12)
         READ(TEXT12,'(3A4)') (HUSE2(I0,NISO2),I0=1,3)
         READ(TEXT12,'(3A4)') (HNAM2(I0,NISO2),I0=1,3)
         CALL XDISET(IPAR,NSTATE,0)
         IPAR(3)=NG
         IPAR(4)=NL
         IPAR(13)=NED+NSPH
         IPAR(19)=NDEL
         IPAR(24)=IDF
      ELSE
*        RECOVER CROSS SECTIONS FROM THE MICROLIB.
         JPEDIT=LCMGID(IPEDIT,'ISOTOPESLIST')
         NISO2=0
         CALL XDISET(ISW,NISOTS,0)
         DO 70 ISO1=1,NISOTS
         IF(IMIX1(ISO1).EQ.IMIL) THEN
            IF(NISOP.GT.0) THEN
               WRITE(TEXT8,'(2A4)') (HUSE1(I0,ISO1),I0=1,2)
               DO 10 JSO=1,NISOP
               IF(NOMISP(JSO).EQ.TEXT8) GO TO 20
   10          CONTINUE
               ISO3=ISO3+1
               ISW(ISO1)=-ISO3
               GO TO 70
            ENDIF
   20       NISO2=NISO2+1
            IF(NISO2.GT.NISOTS) CALL XABORT('COMMIC: NISOTS OVERFLOW.')
            ISW(ISO1)=NISO2
            DO 30 I0=1,2
            HUSE2(I0,NISO2)=HUSE1(I0,ISO1)
   30       CONTINUE
            HUSE2(3,NISO2)=ITEXT
            IF(LISO) HUSE2(3,NISO2)=HUSE1(3,ISO1)
            DO 40 I0=1,3
            HNAM2(I0,NISO2)=HNAM1(I0,ISO1)
   40       CONTINUE
            DENS2(NISO2)=DENS1(ISO1)
            ITYP2(NISO2)=ITYP1(ISO1)
            ITOD2(NISO2)=ITOD1(ISO1)
            VOL2(NISO2)=VOL1(ISO1)
            TEMP2(NISO2)=TEMP1(ISO1)
            KPEDIT=LCMGIL(JPEDIT,ISO1) ! set ISO1-th isotope
            OPCPO=LCMDIL(NPCPO,NISO2) ! set NISO2-th isotope
            CALL LCMEQU(KPEDIT,OPCPO)
*
*           FLUX NORMALIZATION:
            DO 60 IW=1,MIN(NW+1,10)
               WRITE(TEXT12,'(3HNWT,I1)') IW-1
               CALL LCMLEN(OPCPO,TEXT12,ILONG,ITYLCM)
               IF(ILONG.GT.0) THEN
                  CALL LCMGET(OPCPO,TEXT12,WORK)
                  DO 50 IG=1,NG
                  WORK(IG)=WORK(IG)*FNORM
   50             CONTINUE
                  CALL LCMPUT(OPCPO,TEXT12,NG,2,WORK)
               ENDIF
   60       CONTINUE
         ENDIF
   70    CONTINUE
      ENDIF
*----
*  CREATE A NEW MACROSCOPIC RESIDUAL ISOTOPE
*----
      ITRES=0
      ALLOCATE(PYRES(NISO2+1))
      IF(ISO3.GT.0) THEN
         NISO2=NISO2+1
         IF(NISO2.GT.NISOTS) CALL XABORT('COMMIC: NISOTS OVERFLOW(2).')
         CALL LCMOP(IPWORK,'*TEMPORARY*',0,1,0)
         CALL COMRES(IMPX,IPWORK,IPEDIT,NISOTS,NISO2,ISW,FNORM,ITRES,
     1   PYRES)
         OPCPO=LCMDIL(NPCPO,NISO2) ! set NISO2-th isotope
         CALL LCMEQU(IPWORK,OPCPO)
         CALL LCMCL(IPWORK,2)
         TEXT12='*MAC*RES'
         READ(TEXT12,'(3A4)') (HUSE2(I0,NISO2),I0=1,3)
         READ(TEXT12,'(3A4)') (HNAM2(I0,NISO2),I0=1,3)
         DENS2(NISO2)=1.0
         ITYP2(NISO2)=ITRES
         ITOD2(NISO2)=1
         VOL2(NISO2)=VOL2(NISO2-1)
         TEMP2(NISO2)=TEMP2(NISO2-1)
      ENDIF
*----
*  COPY DISCONTINUITY FACTOR INFORMATION AND PERFORM NORMALIZATION
*----
      IF(IDF.NE.0) THEN
         CALL LCMSIX(IPEDIT,'MACROLIB',1)
         CALL LCMLEN(IPEDIT,'ADF',ILONG,ITYLCM)
         IF(ILONG.EQ.0) CALL XABORT('COMMIC: MISSING ADF DIRECTORY IN '
     1   //'EDITION OBJECT.')
         CALL LCMSIX(IPEDIT,'ADF',1)
         CALL LCMSIX(MPCPO,'MACROLIB',1)
         CALL LCMSIX(MPCPO,'ADF',1)
         CALL LCMEQU(IPEDIT,MPCPO)
         IF(IDF.EQ.1)THEN
           ALLOCATE(ADF2(NG,2))
           CALL LCMGET(MPCPO,'ALBS00',ADF2)
           DO IG=1,NG
             ADF2(IG,:2)=ADF2(IG,:2)*FNORM
           ENDDO
           CALL LCMPUT(MPCPO,'ALBS00',NG*2,2,ADF2)
           DEALLOCATE(ADF2)
         ELSEIF(IDF.EQ.2)THEN
           IF(IMPX.GT.5) CALL LCMLIB(MPCPO)
           CALL LCMLEN(MPCPO,'HADF',NTYPE,ITYLCM)
           NTYPE=NTYPE/2
           IF(NTYPE.GT.0) THEN
             ALLOCATE(ADF(NMIL,NG),HADF(NTYPE))
             CALL LCMGTC(MPCPO,'HADF',8,NTYPE,HADF)
             DO ITYPE=1,NTYPE
               CALL LCMLEN(MPCPO,HADF(ITYPE),ILONG,ITYLCM)
               IF(ILONG.NE.NMIL*NG) CALL XABORT('COMMIC: ADF OVERFLOW.')
               CALL LCMGET(MPCPO,HADF(ITYPE),ADF)
               DO IG=1,NG
                 DO IBM=1,NMIL
                   ADF(IBM,IG)=ADF(IBM,IG)*FNORM
                 ENDDO
               ENDDO
               CALL LCMPUT(MPCPO,HADF(ITYPE),NMIL*NG,2,ADF)
             ENDDO
             DEALLOCATE(HADF,ADF)
           ENDIF
         ENDIF
         CALL LCMSIX(MPCPO,' ',2)
         CALL LCMSIX(MPCPO,' ',2)
         CALL LCMSIX(IPEDIT,' ',2)
         CALL LCMSIX(IPEDIT,' ',2)
      ENDIF
*----
*  RECOVER GROUP FORM FACTOR INFORMATION
*----
      IF(NGFF.NE.0) THEN
         IF(.NOT.C_ASSOCIATED(IPEDI2)) CALL XABORT('COMMIC: MISSING ED'
     1   //'ITION OBJECT WITH GFF INFO.')
         CALL COMGFF(MPCPO,IPEDI2,FNORM,NGFF)
      ENDIF
*----
*  RECOVER PHYSICAL ALBEDO INFORMATION
*----
      IF(NALBP.NE.0) THEN
         CALL LCMSIX(IPEDIT,'MACROLIB',1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
         IF(NG.NE.ISTATE(1)) CALL XABORT('COMMIC: INVALID NUMBER OF EN'
     1   //'ERGY GROUPS IN EDITION OBJECT.')
         IF(NALBP.EQ.-1) THEN
            NALBP=ISTATE(8)
         ELSE IF(NALBP.NE.ISTATE(8)) THEN
            CALL XABORT('COMMIC: INVALID NUMBER OF PHYSICAL ALBEDOS IN'
     1      //' EDITION OBJECT.')
         ENDIF
         IF(NALBP.NE.0) THEN
           CALL LCMLEN(IPEDIT,'ALBEDO',ILONG,ITYLCM)
           IF(ILONG.EQ.NALBP*NG) THEN
*            diagonal physical albedos
             ALLOCATE(ALBP(NALBP,NG))
             CALL LCMGET(IPEDIT,'ALBEDO',ALBP)
             CALL LCMSIX(MPCPO,'MACROLIB',1)
             CALL LCMPUT(MPCPO,'ALBEDO',NALBP*NG,2,ALBP)
             CALL LCMSIX(MPCPO,' ',2)
             DEALLOCATE(ALBP)
           ELSE IF(ILONG.EQ.NALBP*NG*NG) THEN
*            matrix physical albedos
             ALLOCATE(ALBP2(NALBP,NG,NG))
             CALL LCMGET(IPEDIT,'ALBEDO',ALBP2)
             CALL LCMSIX(MPCPO,'MACROLIB',1)
             CALL LCMPUT(MPCPO,'ALBEDO',NALBP*NG*NG,2,ALBP2)
             CALL LCMSIX(MPCPO,' ',2)
             DEALLOCATE(ALBP2)
           ELSE
             CALL XABORT('COMMIC: INCONSISTENT ALBEDO INFORMATION.')
           ENDIF
         ENDIF
         CALL LCMSIX(IPEDIT,' ',2)
      ENDIF
*----
*  RESET INFORMATION IN LAMBDA-D, PIFI AND PYIELD
*----
      NDFI2=0
      DO 120 ISO=1,NISO2
      IF(LMACRO.AND.(NDEL.GT.0).AND.(ITYP2(ISO).EQ.2)) THEN
         OPCPO=LCMGIL(NPCPO,ISO) ! set ISO-th isotope
         CALL LCMPUT (OPCPO,'LAMBDA-D',NDEL,2,ZLAMB)
      ELSE IF(ITYP2(ISO).EQ.3) THEN
         OPCPO=LCMGIL(NPCPO,ISO) ! set ISO-th isotope
         CALL LCMLEN(OPCPO,'PIFI',NDFI,ITYLCM)
         IF(NDFI.GT.0) THEN
           ALLOCATE(IPIFI(NDFI),PYIELD(NDFI),IPIFI2(NDFI+1),
     1     PYIEL2(NDFI+1))
           CALL LCMGET(OPCPO,'PIFI',IPIFI)
           CALL LCMGET(OPCPO,'PYIELD',PYIELD)
           NDFI2=0
           LRES=.FALSE.
           DO 110 I=1,NDFI
             IFI=IPIFI(I)
             IF(IFI.GT.NISOTS) CALL XABORT('COMMIC: NISOTS OVERFLOW.')
             IF(ISW(IFI).GT.0) THEN
                NDFI2=NDFI2+1
                IPIFI2(NDFI2)=ISW(IFI)
                PYIEL2(NDFI2)=PYIELD(I)
             ELSE IF(ISW(IFI).LT.0) THEN
                LRES=.TRUE.
             ENDIF
  110      ENDDO
           IF(LRES) THEN
             NDFI2=NDFI2+1
             IPIFI2(NDFI2)=NISO2
             PYIEL2(NDFI2)=PYRES(ISO)
           ENDIF
           IF(NDFI2.GT.0) THEN
             CALL LCMPUT(OPCPO,'PIFI',NDFI2,1,IPIFI2)
             CALL LCMPUT(OPCPO,'PYIELD',NDFI2,2,PYIEL2)
           ENDIF
           DEALLOCATE(PYIEL2,IPIFI2,PYIELD,IPIFI)
         ENDIF
      ENDIF
  120 CONTINUE
      DEALLOCATE(PYRES)
*
      IPAR(1)=1
      IPAR(2)=NISO2
      IPAR(11)=0
      IPAR(14)=1
      IPAR(20)=NDFI2
      IPAR(25)=NW
      TEXT12='L_LIBRARY'
      CALL LCMPTC(MPCPO,'SIGNATURE',12,1,TEXT12)
      CALL LCMPUT(MPCPO,'STATE-VECTOR',NSTATE,1,IPAR)
      CALL XDISET(ISW,NISO2,1)
      CALL LCMPUT(MPCPO,'ISOTOPESMIX',NISO2,1,ISW)
      CALL LCMPUT(MPCPO,'ISOTOPESUSED',3*NISO2,3,HUSE2)
      CALL LCMPUT(MPCPO,'ISOTOPERNAME',3*NISO2,3,HNAM2)
      CALL LCMPUT(MPCPO,'ISOTOPESDENS',NISO2,2,DENS2)
      CALL LCMPUT(MPCPO,'ISOTOPESTYPE',NISO2,1,ITYP2)
      CALL LCMPUT(MPCPO,'ISOTOPESTODO',NISO2,1,ITOD2)
      CALL LCMPUT(MPCPO,'ISOTOPESVOL',NISO2,2,VOL2)
      CALL LCMPUT(MPCPO,'ISOTOPESTEMP',NISO2,2,TEMP2)
      CALL LCMPUT(MPCPO,'MIXTURESVOL',1,2,VOLMIX(IMIL))
      CALL LCMPUT(MPCPO,'K-EFFECTIVE',1,2,EIGENK)
      CALL LCMPUT(MPCPO,'K-INFINITY',1,2,EIGINF)
      IF(B2.NE.0.0) CALL LCMPUT(MPCPO,'B2  B1HOM',1,2,B2)
      IF((NED+NSPH).GT.0) CALL LCMPUT(MPCPO,'ADDXSNAME-P0',2*(NED+NSPH),
     1 3,HVECT)
      CALL LCMPUT(MPCPO,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(MPCPO,'DELTAU',NG,2,DELT)
      IF(IMPX.GT.2) WRITE(6,140) IMIL,NISO2
      IF(IMPX.GT.5) CALL LCMLIB(MPCPO)
      IF(LMACRO) DEALLOCATE(ZLAMB)
  130 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOLMIX,DELT,ENER,WORK,VOL2,TEMP2,DENS2,VOL1,TEMP1,
     1 DENS1)
      DEALLOCATE(HVECT,ISW,ITOD2,ITYP2,HNAM2,HUSE2,ITOD1,ITYP1,IMIX1,
     1 HNAM1,HUSE1)
      RETURN
*
  140 FORMAT(39H COMMIC: PROCESSING HOMOGENIZED MIXTURE,I4,9H CONTAINI,
     1 2HNG,I5,10H ISOTOPES.)
      END
