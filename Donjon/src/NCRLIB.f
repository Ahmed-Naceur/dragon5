*DECK NCRLIB
      SUBROUTINE NCRLIB(MAXNIS,MAXISO,MAXFEL,IPLIB,IPCPO,IACCS,NMIL,
     1 NMIX,NGRP,NGFF,NALBP,IDF,IMPX,NCAL,TERP,NISO,LISO,HISO,CONC,
     2 MIXC,LXS,LRES,LPURE,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the microlib by scanning the NCAL elementary calculations and
* weighting them with TERP factors.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* MAXNIS  maximum value of NISO(I) in user data.
* MAXISO  maximum allocated space for output microlib TOC information.
* MAXFEL  number of fuel rings used for the micro-depletion.
* IPLIB   address of the output microlib LCM object.
* IPCPO   address of the multicompo object.
* IACCS   =0 microlib is created; =1 ... is updated.
* NMIL    number of material mixtures in the multicompo.
* NMIX    maximum number of material mixtures in the microlib.
* NGRP    number of energy groups.
* NGFF    number of group form factors per energy group.
* NALBP   number of physical albedos per energy group.
* IDF     ADF type, 0 = none, 1 = Albedo, 2 = FD_B/FD_C/...
* IMPX    print parameter (equal to zero for no print).
* NCAL    number of elementary calculations in the multicompo.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes.
*         A value of -99.99 is set to indicate that the multicompo value
*         is used.
* MIXC    mixture index in the multicompo corresponding to each microlib
*         mixture. Equal to zero if a microlib mixture is not updated.
* LXS     =.true. if keyword 'ALLX' is specified
* LRES    =.true. if the interpolation is done without updating isotopic 
*         densities
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
* B2      buckling
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPCPO
      INTEGER MAXNIS,MAXISO,MAXFEL,IACCS,NMIL,NMIX,NGRP,NGFF,NALBP,
     1 IDF,IMPX,NCAL,NISO(NMIX),HISO(2,NMIX,MAXNIS),MIXC(NMIX)
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXNIS),B2
      LOGICAL LISO(NMIX),LXS,LRES,LPURE
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXED=50
      INTEGER, PARAMETER::NSTATE=40
      INTEGER I0, IACCOLD, IBMOLD, IBM, ICAL, IED1, IED2, IGR,  ILONG,
     & ISO, ITRANC, ITYLCM, I, JSO1, JSO, J, KSO1, KSO, NBISO1, NBISO2,
     & NBISOT2, NBISOT, NBRG, NCOMB2, NCOMB, NDEL, NDEPL, NDFI, NED1,
     & NED2, NFINF, NL, NW, NTYPE
      REAL WEIGHT
      CHARACTER TEXT12*12,HNAME*12,HSMG*131,HVECT1(MAXED)*8,
     1 HVECT2(MAXED)*8,CHAR1*4,CHAR2*4,HHISO*8
      INTEGER ISTATE(NSTATE)
      LOGICAL LUSER
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO,MPCPO,JPLIB,KPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYP1,ITOD1,IMIX2,ITYP2,
     1 ITOD2,MILVO,IMICR
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HUSE1,HNAM1,HUSE2,HNAM2
      REAL, ALLOCATABLE, DIMENSION(:) :: TEMP1,VOL1,DENS2,TEMP2,VOL2,
     1 DENS3,TEMP3,VOL3,ENER,DELT,VOLMI2,GAR1,GAR2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FACT,DENS1
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPLIST
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(HUSE1(3,MAXISO),HNAM1(3,MAXISO),ITYP1(MAXISO),
     1 ITOD1(MAXISO),IMIX2(MAXISO),ITYP2(MAXISO),ITOD2(MAXISO),
     2 HUSE2(3,MAXISO),HNAM2(3,MAXISO),MILVO(NMIX))
      ALLOCATE(TEMP1(MAXISO),VOL1(MAXISO),DENS2(MAXISO),TEMP2(MAXISO),
     1 VOL2(MAXISO),ENER(NGRP+1),DELT(NGRP),VOLMI2(NMIX),IPLIST(MAXISO))
      IACCOLD=IACCS ! for ADF and GFF 
*----
*  MICROLIB INITIALIZATION
*----
      ITRANC=0
      VOLMI2(:NMIX)=0.0
      DENS2(:MAXISO)=0.0
      VOL2(:MAXISO)=0.0
      TEMP2(:MAXISO)=0.0
      IMIX2(:MAXISO)=0
      ITYP2(:MAXISO)=0
      ITOD2(:MAXISO)=0
      IPLIST(:MAXISO)=C_NULL_PTR
      IF(IACCS.EQ.0) THEN
         IF(LRES) CALL XABORT('NCRLIB: RES OPTION IS INVALID.')
         NBISO2=0
         NCOMB2=0
         NED2=0
         TEXT12='L_LIBRARY'
         CALL LCMPTC(IPLIB,'SIGNATURE',12,1,TEXT12)
      ELSE
         CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NMIX) CALL XABORT('NCRLIB: INVALID NUMBER OF '
     1   //'MATERIAL MIXTURES IN THE MICROLIB.')
         IF(ISTATE(3).NE.NGRP) CALL XABORT('NCRLIB: INVALID NUMBER OF '
     1   //'ENERGY GROUPS IN THE MICROLIB.')
         NBISO2=ISTATE(2)
         NCOMB2=ISTATE(12)
         IF(NBISO2.GT.MAXISO) CALL XABORT('NCRLIB: MAXISO OVERFLOW(1).')
         NED2=ISTATE(13)
         IF(NED2.GT.MAXED) CALL XABORT('NCRLIB: MAXED OVERFLOW.')
         CALL LCMLEN(IPLIB,'MIXTURESVOL',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
           CALL LCMGET(IPLIB,'MIXTURESVOL',VOLMI2)
         ELSE
           VOLMI2(:NMIX)=0.0
         ENDIF
         CALL LCMGET(IPLIB,'ISOTOPESUSED',HUSE2)
         CALL LCMGET(IPLIB,'ISOTOPERNAME',HNAM2)
         CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYP2)
         CALL LCMGET(IPLIB,'ISOTOPESTODO',ITOD2)
         CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS2)
         CALL LCMGET(IPLIB,'ISOTOPESVOL',VOL2)
         CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX2)
         CALL LCMGET(IPLIB,'ISOTOPESTEMP',TEMP2)
         IF(NED2.GT.0) CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
         CALL LCMGET(IPLIB,'ENERGY',ENER)
         CALL LCMGET(IPLIB,'DELTAU',DELT)
      ENDIF
*----
*  RECOVER NDEPL
*----
      NDEPL=0
      CALL LCMLEN(IPCPO,'DEPL-CHAIN',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPCPO,'DEPL-CHAIN',1)
        CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
        NDEPL=ISTATE(1)
        CALL LCMSIX(IPCPO,' ',2)
      ENDIF
*----
*  LOOP OVER MICROLIB MIXTURES
*----
      ALLOCATE(DENS3(MAXISO),TEMP3(MAXISO),VOL3(MAXISO))
      MILVO(:NMIX)=0
      NCOMB=0
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      DO 190 IBM=1,NMIX
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.EQ.0) GO TO 190
      IF(NISO(IBM).GT.MAXNIS) CALL XABORT('NCRLIB: MAXNIS OVERFLOW.')
      KPCPO=LCMGIL(JPCPO,IBMOLD)
      LPCPO=LCMGID(KPCPO,'CALCULATIONS')
*----
*  FIND THE VALUE OF NBISO1 IN MIXTURE IBM
*----
      DO ICAL=1,NCAL
        IF(TERP(ICAL,IBM).EQ.0.0) CYCLE
        MPCPO=LCMGIL(LPCPO,ICAL)
        CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
        NBISO1=ISTATE(2)
        CALL LCMGET(MPCPO,'ISOTOPESUSED',HUSE1)
        CALL LCMGET(MPCPO,'ISOTOPERNAME',HNAM1)
        EXIT
      ENDDO
      ALLOCATE(FACT(NCAL,NBISO1),DENS1(NBISO1,NCAL))
*----
*  LOOP OVER ELEMENTARY CALCULATIONS
*----
      JSO1=0
      DENS3(:NBISO1)=0.0
      VOL3(:NBISO1)=0.0
      TEMP3(:NBISO1)=0.0
      DO 50 ICAL=1,NCAL
      WEIGHT=TERP(ICAL,IBM)
      IF(WEIGHT.EQ.0.0) GO TO 50
      MPCPO=LCMGIL(LPCPO,ICAL)
      IF(IMPX.GT.0) THEN
         WRITE(IOUT,'(38H NCRLIB: MULTICOMPO ACCESS FOR MIXTURE,I8,
     1   5H (<==,I4,17H) AND CALCULATION,I8,9H. WEIGHT=,1P,E12.4)')
     2   IBM,IBMOLD,ICAL,WEIGHT
         IF(IMPX.GT.50) CALL LCMLIB(MPCPO)
      ENDIF
      CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.1) CALL XABORT('NCRLIB: INVALID NUMBER OF MATERI'
     1 //'AL MIXTURES IN THE MULTICOMPO.')
      IF(ISTATE(2).NE.NBISO1) CALL XABORT('NCRLIB: INVALID NBISO1.')
      IF(ISTATE(3).NE.NGRP) CALL XABORT('NCRLIB: INVALID NUMBER OF ENE'
     1 //'RGY GROUPS IN THE MULTICOMPO.')
      NL=ISTATE(4)
      ITRANC=ISTATE(5)
      NDEPL=MAX(ISTATE(11),NDEPL)
      NED1=ISTATE(13)
      NDEL=ISTATE(19)
      NDFI=ISTATE(20)
      NW=ISTATE(25)
      IF(NED1.GT.MAXED) CALL XABORT('NCRLIB: MAXED OVERFLOW.')
      CALL LCMLEN(MPCPO,'MIXTURESVOL',ILONG,ITYLCM)
      IF(ILONG.GT.0) CALL LCMGET(MPCPO,'MIXTURESVOL',VOLMI2(IBM))
      CALL LCMGET(MPCPO,'ISOTOPESDENS',DENS1(1,ICAL))
      CALL LCMGET(MPCPO,'ISOTOPESTYPE',ITYP1)
      CALL LCMGET(MPCPO,'ISOTOPESTODO',ITOD1)
      CALL LCMGET(MPCPO,'ISOTOPESVOL',VOL1)
      CALL LCMGET(MPCPO,'ISOTOPESTEMP',TEMP1)
      IF(NED1.GT.0) CALL LCMGTC(MPCPO,'ADDXSNAME-P0',8,NED1,HVECT1)
      CALL LCMGET(MPCPO,'ENERGY',ENER)
      CALL LCMGET(MPCPO,'DELTAU',DELT)
      DO 30 IED1=1,NED1
      DO 20 IED2=1,NED2
      IF(HVECT1(IED1).EQ.HVECT2(IED2)) GO TO 30
   20 CONTINUE
      NED2=NED2+1
      HVECT2(NED2)=HVECT1(IED1)
   30 CONTINUE
      DO 49 ISO=1,NBISO1 ! multicompo isotope
      WRITE(TEXT12,'(2A4)') (HUSE1(I,ISO),I=1,2)
      IF(LXS) THEN
        WRITE(CHAR1,'(A4)') HUSE1(3,ISO)
        READ(CHAR1,'(I4)') NBRG
        NBISOT=NBRG+MAXFEL*(IBM-1)
        IF(NBISOT.GT.9999) CALL XABORT('NCRLIB: NBISOT OVERFLOW.')
        WRITE(TEXT12,'(2A4,I4.4)') (HUSE1(I,ISO),I=1,2),NBISOT
      ENDIF
      KSO1=0
      DO 40 KSO=1,NISO(IBM) ! user-selected isotope
      WRITE(HHISO,'(2A4)') (HISO(I0,IBM,KSO),I0=1,2)
      IF(TEXT12(:8).EQ.HHISO) THEN
        KSO1=KSO
        GO TO 45
      ENDIF
   40 CONTINUE
   45 LUSER=.FALSE.
      IF(KSO1.GT.0) LUSER=(CONC(IBM,KSO1).NE.-99.99)
      IF(LUSER) DENS1(ISO,ICAL)=CONC(IBM,KSO1)
      DENS3(ISO)=DENS3(ISO)+WEIGHT*DENS1(ISO,ICAL)
      VOL3(ISO)=VOL3(ISO)+WEIGHT*VOL1(ISO)
      TEMP3(ISO)=TEMP3(ISO)+WEIGHT*TEMP1(ISO)
   49 CONTINUE
   50 CONTINUE
      FACT(:NCAL,:NBISO1)=1.0
      IF(.NOT.LPURE) THEN
        DO ICAL=1,NCAL
          IF(TERP(ICAL,IBM).EQ.0.0) CYCLE
          DO ISO=1,NBISO1 ! multicompo isotope
            IF(DENS3(ISO).GT.DENS1(ISO,ICAL)*1.0E-9) THEN
              FACT(ICAL,ISO)=DENS1(ISO,ICAL)/DENS3(ISO)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(DENS1)
*----
*  SET NUMBER DENSITIES AND VOLUMES IN OUTPUT MICROLIB
*----
      IF(LRES) THEN
*       -- Number densities are left unchanged except if they are
*       -- listed in HISO array.
        DO 60 KSO=1,NISO(IBM) ! user-selected isotope
          DO JSO=1,NBISO2 ! microlib isotope
            IF(IMIX2(JSO).NE.IBM) CYCLE
            IF((HISO(1,IBM,KSO).EQ.HUSE2(1,JSO)).AND.
     1         (HISO(2,IBM,KSO).EQ.HUSE2(2,JSO))) THEN
              IF(CONC(IBM,KSO).EQ.-99.99) THEN
*               -- Only number densities of isotopes set with "MICR" and
*               -- "*" keywords are interpolated
                DENS2(JSO)=0.0
                DO ISO=1,NBISO1 ! multicompo isotope
                  JSO1=0
                  IF((HUSE1(1,ISO).EQ.HUSE2(1,JSO)).AND.
     1               (HUSE1(2,ISO).EQ.HUSE2(2,JSO))) THEN
                    IF(ITYP1(ISO).NE.ITYP2(JSO)) THEN
                      WRITE(HSMG,500) 'ITYP',ISO,ITYP1(ISO),ITYP2(JSO)
                      CALL XABORT(HSMG)
                    ENDIF
                    IF(LXS) THEN
                      WRITE(CHAR1,'(A4)') HUSE1(3,ISO)
                      WRITE(CHAR2,'(A4)') HUSE2(3,JSO)
                      READ(CHAR2,'(I4.4)') NBISOT2
                      NBISOT2=NBISOT2-MAXFEL*(IBM-1)
                      WRITE(CHAR2,'(I4.4)') NBISOT2
                      IF(CHAR1.EQ.CHAR2) THEN
                        JSO1=JSO
                        GO TO 55
                      ENDIF
                    ELSE
                      JSO1=JSO
                      GO TO 55
                    ENDIF
   55               IF(JSO1.EQ.0) CALL XABORT('NCRLIB: JSO1=0')
                    DENS2(JSO1)=DENS2(JSO1)+DENS3(ISO)
                    TEMP2(JSO1)=TEMP3(ISO)
                  ENDIF
                ENDDO
              ELSE IF(CONC(IBM,KSO).NE.-99.99) THEN
*               -- Number densities of isotopes set with "MICR" and
*               -- fixed value are forced to this value
                DENS2(JSO)=CONC(IBM,KSO)
              ENDIF
              GO TO 60
            ENDIF
          ENDDO
          WRITE(HSMG,'(31HNCRLIB: UNABLE TO FIND ISOTOPE ,2A4,6H IN MI,
     1    5HXTURE,I8,1H.)') HISO(1,IBM,KSO),HISO(2,IBM,KSO),IBM
          CALL XABORT(HSMG)
   60   CONTINUE
      ELSE
*       -- Number densities are interpolated or not according to
*       -- ALL/ONLY option
        DO JSO=1,NBISO2 ! microlib isotope
          IF(IBM.EQ.IMIX2(JSO)) THEN
            DO ISO=1,NBISO1 ! multicompo isotope
              IF((HUSE1(1,ISO).EQ.HUSE2(1,JSO)).AND.
     1           (HUSE1(2,ISO).EQ.HUSE2(2,JSO))) THEN
                DENS2(JSO)=0.0
                VOL2(JSO)=0.0
                CYCLE
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO 110 ISO=1,NBISO1 ! multicompo isotope
        WRITE(TEXT12,'(2A4)') (HUSE1(I,ISO),I=1,2)
        IF(LXS) THEN
          WRITE(CHAR1,'(A4)') HUSE1(3,ISO)
          READ(CHAR1,'(I4)') NBRG
          NBISOT=NBRG+MAXFEL*(IBM-1)
          IF(NBISOT.GT.9999) CALL XABORT('NCRLIB: NBISOT OVERFLOW.')
          WRITE(TEXT12,'(2A4,I4.4)') (HUSE1(I,ISO),I=1,2),NBISOT
        ENDIF
        IF(.NOT.LISO(IBM)) THEN
*         --ONLY option
          DO KSO=1,NISO(IBM) ! user-selected isotope
            WRITE(HHISO,'(2A4)') (HISO(I0,IBM,KSO),I0=1,2)
            IF(TEXT12(:8).EQ.HHISO) GO TO 65
          ENDDO
          GO TO 110
        ENDIF
   65   DO 70 JSO=1,NBISO2 ! microlib isotope
        JSO1=0
        IF((HUSE1(1,ISO).EQ.HUSE2(1,JSO)).AND.(HUSE1(2,ISO).EQ.
     1  HUSE2(2,JSO)).AND.(IMIX2(JSO).EQ.IBM)) THEN
          IF(ITYP1(ISO).NE.ITYP2(JSO)) THEN
            WRITE(HSMG,500) 'ITYP',ISO,ITYP1(ISO),ITYP2(JSO)
            CALL XABORT(HSMG)
          ENDIF
          IF(LXS) THEN
            WRITE(CHAR1,'(A4)') HUSE1(3,ISO)
            WRITE(CHAR2,'(A4)') HUSE2(3,JSO)
            READ(CHAR2,'(I4.4)') NBISOT2
            NBISOT2=NBISOT2-MAXFEL*(IBM-1)
            WRITE(CHAR2,'(I4.4)') NBISOT2
            IF(CHAR1.EQ.CHAR2) THEN
              JSO1=JSO
              GO TO 100
            ENDIF
          ELSE
            JSO1=JSO
            GO TO 100
          ENDIF
        ENDIF
   70   CONTINUE
        NBISO2=NBISO2+1
        IF(NBISO2.GT.MAXISO) THEN
          WRITE(IOUT,'(/16H NCRLIB: NBISO2=,I6,8H MAXISO=,I6)') NBISO2,
     1    MAXISO
          CALL XABORT('NCRLIB: MAXISO OVERFLOW(2).')
        ENDIF
        READ(TEXT12,'(3A4)') (HUSE2(I0,NBISO2),I0=1,3)
        DO 80 I0=1,3
        HNAM2(I0,NBISO2)=HNAM1(I0,ISO)
   80   CONTINUE
        IMIX2(NBISO2)=IBM
        ITYP2(NBISO2)=ITYP1(ISO)
        ITOD2(NBISO2)=ITOD1(ISO)
        IF(ITYP2(NBISO2).EQ.1) ITOD2(NBISO2)=1
        JSO1=NBISO2
        IF(ITOD2(NBISO2).NE.1) THEN
          DO 90 J=1,NCOMB
          IF(IBM.EQ.MILVO(J)) GO TO 100
   90     CONTINUE
          NCOMB=NCOMB+1
          IF(NCOMB.GT.NMIX) CALL XABORT('NCRLIB: MILVO OVERFLOW.')
          MILVO(NCOMB)=IBM
        ENDIF
  100   DENS2(JSO1)=DENS2(JSO1)+DENS3(ISO)
        VOL2(JSO1)=VOL2(JSO1)+VOL3(ISO)
        TEMP2(JSO1)=TEMP3(ISO)
  110   CONTINUE
      ENDIF
*----
*  SELECT MICROLIB ISOTOPES CORRESPONDING TO MULTICOMPO ISOTOPES
*----
      ALLOCATE(IMICR(NBISO1))
      IMICR(:NBISO1)=0
      DO 130 ISO=1,NBISO2 ! microlib isotope
      IF(IMIX2(ISO).NE.IBM) GO TO 130
      DO 120 JSO=1,NBISO1 ! multicompo isotope
      IF((HUSE1(1,JSO).EQ.HUSE2(1,ISO)).AND.(HUSE1(2,JSO).EQ.
     1 HUSE2(2,ISO))) THEN
        IF(LXS) THEN
          WRITE(CHAR1,'(A4)') HUSE1(3,JSO)
          WRITE(CHAR2,'(A4)') HUSE2(3,ISO)
          READ(CHAR1,'(I4.4)') NBRG
          NBISOT=NBRG+MAXFEL*(IBM-1)
          READ(CHAR2,'(I4.4)') NBISOT2
          IF(NBISOT.EQ.NBISOT2) THEN
            IMICR(JSO)=ISO
            GO TO 130
          ENDIF
        ELSE
          IMICR(JSO)=ISO
          GO TO 130
        ENDIF
      ENDIF
  120 CONTINUE
      WRITE(TEXT12,'(3A4)') (HUSE2(I0,ISO),I0=1,3)
      CALL XABORT('NCRLIB: UNABLE TO FIND '//TEXT12//'.')
  130 CONTINUE
*----
*  PROCESS ISOTOPE DIRECTORIES FOR MICROLIB MIXTURE IBM
*----
      DO 180 JSO=1,NBISO1 ! multicompo isotope
      ISO=IMICR(JSO) ! microlib isotope
      IF(ISO.EQ.0) GO TO 180
      WRITE(HNAME,'(3A4)') (HUSE1(I0,JSO),I0=1,3)
      CALL LCMOP(KPLIB,'*ISOTOPE*',0,1,0)
      IPLIST(ISO)=KPLIB ! set isot ISO
      CALL NCRISO(KPLIB,LPCPO,NBISO1,IMICR,HNAME,JSO,IBMOLD,NCAL,NGRP,
     1 NL,NW,NED2,HVECT2,NDEL,NDFI,IMPX,FACT(1,JSO),TERP(1,IBM),LPURE)
  180 CONTINUE
      DEALLOCATE(IMICR,FACT)
  190 CONTINUE
      DEALLOCATE(VOL3,TEMP3,DENS3)
      DEALLOCATE(VOL1,TEMP1,MILVO,ITOD1,ITYP1,HNAM1,HUSE1)
*----
*  CREATE ISOTOPE LIST DIRECTORY IN MICROLIB
*----
      JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',NBISO2)
      DO 195 ISO=1,NBISO2 ! microlib isotope
      IF(C_ASSOCIATED(IPLIST(ISO))) THEN
        KPLIB=LCMDIL(JPLIB,ISO) ! step up isot ISO
        CALL LCMEQU(IPLIST(ISO),KPLIB)
        CALL LCMCL(IPLIST(ISO),2)
      ENDIF
  195 CONTINUE
      DEALLOCATE(IPLIST)
*----
*  MICROLIB FINALIZATION
*----
      IF(.NOT.LRES) THEN
        ISTATE(:NSTATE)=0
        ISTATE(1)=NMIX
        ISTATE(2)=NBISO2
        ISTATE(3)=NGRP
        ISTATE(4)=NL
        ISTATE(5)=ITRANC
        ISTATE(7)=1
        ISTATE(11)=NDEPL
        ISTATE(12)=NCOMB+NCOMB2
        ISTATE(13)=NED2
        ISTATE(14)=NMIX
        ISTATE(18)=1
        ISTATE(19)=NDEL
        ISTATE(20)=NDFI
        ISTATE(22)=MAXISO/NMIX
        IF(NBISO2.EQ.0) CALL XABORT('NCRLIB: NBISO2=0.')
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        CALL LCMPUT(IPLIB,'MIXTURESVOL',NMIX,2,VOLMI2)
        CALL LCMPUT(IPLIB,'ISOTOPESUSED',3*NBISO2,3,HUSE2)
        CALL LCMPUT(IPLIB,'ISOTOPERNAME',3*NBISO2,3,HNAM2)
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESMIX',NBISO2,1,IMIX2)
        CALL LCMPUT(IPLIB,'ISOTOPESTYPE',NBISO2,1,ITYP2)
        CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO2,1,ITOD2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
        CALL LCMPUT(IPLIB,'ISOTOPESTEMP',NBISO2,2,TEMP2)
        IF(NED2.GT.0) CALL LCMPTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
        CALL LCMPUT(IPLIB,'ENERGY',NGRP+1,2,ENER)
        CALL LCMPUT(IPLIB,'DELTAU',NGRP,2,DELT)
      ELSE
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
        CALL LCMPUT(IPLIB,'ISOTOPESTEMP',NBISO2,2,TEMP2)
      ENDIF
      IF(IMPX.GT.5) CALL LCMLIB(IPLIB)
      IACCS=1
      DEALLOCATE(VOLMI2,DELT,ENER,VOL2,TEMP2,DENS2,HNAM2,HUSE2,ITOD2,
     1 ITYP2,IMIX2)
*----
*  BUILD EMBEDDED MACROLIB
*----
      CALL SPHEMB(IPLIB,IPCPO,NGRP,NMIX,MIXC)
*----
*  INCLUDE LEAKAGE IN THE MACROLIB (USED ONLY FOR NON-REGRESSION TESTS)
*----
      IF(B2.NE.0.0) THEN
        IF(IMPX.GT.0) WRITE(6,'(/34H NCRLIB: INCLUDE LEAKAGE IN THE MA,
     1  11HCROLIB (B2=,1P,E12.5,2H).)') B2
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        JPLIB=LCMGID(IPLIB,'GROUP')
        ALLOCATE(GAR1(NMIX),GAR2(NMIX))
        DO 210 IGR=1,NGRP
          KPLIB=LCMGIL(JPLIB,IGR)
          CALL LCMGET(KPLIB,'NTOT0',GAR1)
          CALL LCMGET(KPLIB,'DIFF',GAR2)
          DO 200 IBM=1,NMIX
            IF(MIXC(IBM).NE.0) GAR1(IBM)=GAR1(IBM)+B2*GAR2(IBM)
  200     CONTINUE
          CALL LCMPUT(KPLIB,'NTOT0',NMIX,2,GAR1)
  210   CONTINUE
        DEALLOCATE(GAR2,GAR1)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  PROCESS ADF, GFF and physical albedos (if required)
*----
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      KPCPO=LCMGIL(JPCPO,1)
      LPCPO=LCMGID(KPCPO,'CALCULATIONS')
      MPCPO=LCMGIL(LPCPO,1)
      CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
      IDF=ISTATE(24)
      NTYPE=0
      IF(IDF.EQ.1) THEN
        NTYPE=2
      ELSE IF(IDF.EQ.2) THEN
        CALL LCMSIX(MPCPO,'MACROLIB',1)
        CALL LCMLEN(MPCPO,'ADF',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CALL XABORT('NCRLIB: MISSING ADF DIRECTORY I'
     1  //'N MULTICOMPO OBJECT.')
        CALL LCMSIX(MPCPO,'ADF',1)
        CALL LCMGET(MPCPO,'NTYPE',NTYPE)
        CALL LCMSIX(MPCPO,' ',2)
        CALL LCMSIX(MPCPO,' ',2)
      ENDIF
      IF(NGFF.GT.0) THEN
        CALL LCMSIX(MPCPO,'MACROLIB',1)
        CALL LCMLEN(MPCPO,'GFF',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CALL XABORT('NCRLIB: MISSING GFF DIRECTORY I'
     1  //'N MULTICOMPO OBJECT.')
        CALL LCMSIX(MPCPO,'GFF',1)
        CALL LCMLEN(MPCPO,'FINF_NUMBER ',NFINF,ITYLCM)
        CALL LCMSIX(MPCPO,' ',2)
        CALL LCMSIX(MPCPO,' ',2)
      ENDIF
      IF((NGFF.GT.0).OR.(IDF.GT.0).OR.(NALBP.GT.0))THEN
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL NCRAGF(IPLIB,IPCPO,IACCOLD,NMIL,NMIX,NGRP,NGFF,NALBP,IMPX,
     1  NCAL,TERP,MIXC,IDF,NTYPE,NFINF)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
      RETURN
*
  500 FORMAT(8HNCRLIB: ,A,1H(,I4,2H)=,2I5)
      END
