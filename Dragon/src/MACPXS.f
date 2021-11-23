*DECK MACPXS
      SUBROUTINE MACPXS(IPLIST,MAXFIS,NGROUP,NBMIX,NIFISS,NANISO,NDELG,
     >                  ITRANC,LNEWXS,XSTOTL,XSTOT1,XSFISS,XSSPEC,
     >                  XSFIXE,XSTRAN,XSDIFF,XSNFTO,XSH,XSSCAT,NEDMAC,
     >                  ISCATA,XSNUDL,XSCHDL,XSDIFX,XSDIFY,XSDIFZ,
     >                  XSOVRV,XSINT0,XSINT1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transfer cross section information on the macrolib.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIST  LCM pointer to the macrolib.
* MAXFIS  set to max(1,NIFISS).
* NGROUP  number of energy groups.
* NBMIX   maximum number of mixtures.
* NIFISS  number of fissile isotopes.
* NANISO  maximum Legendre order (=1 isotropic collision;
*         =2 linearly anisotropic collision).
* NDELG   number of precursor groups for delayed neutrons.
* ITRANC  transport correction option (=0 no correction; =1 Apollo-
*         type; =2 recover TRANC record; =4 leakage correction alone).
* LNEWXS  check change in cross sections.
* XSTOTL  P0 total cross section of mixture.
* XSTOT1  P1 total cross section of mixture.
* XSFISS  nu*fission cross section of mixture.
* XSSPEC  fission spectrum.
* XSFIXE  fixe sources.
* XSTRAN  transport correction.
* XSDIFF  isotropic diffusion coefficient.
* XSNFTO  fission cross section of mixture.
* XSH     power factor.
* XSSCAT  scattering cross section of mixture/group.
* NEDMAC  number of macro edit cross sections.
* XSNUDL  delayed nu*fission cross section of mixture.
* XSCHDL  delayed-neutron fission spectrum.
* XSDIFX  x-directed diffusion coefficients.
* XSDIFY  y-directed diffusion coefficients.
* XSDIFZ  z-directed diffusion coefficients.
* XSOVRV  reciprocal neutron velocities.
* XSINT0  P0 volume-integrated flux of mixture.
* XSINT1  P1 volume-integrated flux of mixture.
* ISCATA  check for scattering anisotropy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER   MAXFIS,NGROUP,NBMIX,NIFISS,NANISO,NDELG,ITRANC,NEDMAC,
     >          ISCATA(NANISO)
      REAL      XSTOTL(NBMIX,NGROUP),XSTOT1(NBMIX,NGROUP),
     >          XSFISS(NBMIX,MAXFIS,NGROUP),XSSPEC(NBMIX,MAXFIS,NGROUP),
     >          XSFIXE(NBMIX,NGROUP),XSTRAN(NBMIX,NGROUP),
     >          XSDIFF(NBMIX,NGROUP),XSNFTO(NBMIX,NGROUP),
     >          XSH(NBMIX,NGROUP),XSSCAT(NGROUP,NBMIX,NANISO,NGROUP),
     >          XSNUDL(NBMIX,MAXFIS,NDELG,NGROUP),
     >          XSCHDL(NBMIX,MAXFIS,NDELG,NGROUP),
     >          XSDIFX(NBMIX,NGROUP),XSDIFY(NBMIX,NGROUP),
     >          XSDIFZ(NBMIX,NGROUP),XSOVRV(NBMIX,NGROUP),
     >          XSINT0(NBMIX,NGROUP),XSINT1(NBMIX,NGROUP)
      LOGICAL   LNEWXS(18)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXNED=50)
      TYPE(C_PTR) JPLIST,KPLIST
      CHARACTER CANISO*2,CHID*12,NUSIGD*12,HVECT(MAXNED)*8
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGSCT,IFGSCT,IPOSCT
      REAL, ALLOCATABLE, DIMENSION(:) :: XSWORK
*----
*  SCRATCH STORAGE ALLOCATION
*   INGSCT  number of scattering group for cross sections.
*   IFGSCT  first scattering group for cross sections.
*   IPOSCT  material position in scattering.
*----
      ALLOCATE(INGSCT(NBMIX),IFGSCT(NBMIX),IPOSCT(NBMIX))
      ALLOCATE(XSWORK(NBMIX*(2+NGROUP)))
*----
*  GET NUMBER OF MACRO EDIT X-SECTIONS
*----
      CALL LCMLEN(IPLIST,'ADDXSNAME-P0',ILCMLN,ILCMTY)
      NEDMAC=ILCMLN/2
      IF(NEDMAC.GT.MAXNED) CALL XABORT('MACPXS: MAXNED OVERFLOW(1).')
      IF(NEDMAC.GT.0) CALL LCMGTC(IPLIST,'ADDXSNAME-P0',8,NEDMAC,HVECT)
      IF(LNEWXS(8)) THEN
        DO IED=1,NEDMAC
          IF(HVECT(IED).EQ.'H-FACTOR') GO TO 5
        ENDDO
        NEDMAC=NEDMAC+1
        IF(NEDMAC.GT.MAXNED) CALL XABORT('MACPXS: MAXNED OVERFLOW(2).')
        HVECT(NEDMAC)='H-FACTOR'
      ENDIF
    5 IF(LNEWXS(15)) THEN
        DO IED=1,NEDMAC
          IF(HVECT(IED).EQ.'OVERV') GO TO 10
        ENDDO
        NEDMAC=NEDMAC+1
        IF(NEDMAC.GT.MAXNED) CALL XABORT('MACPXS: MAXNED OVERFLOW(3).')
        HVECT(NEDMAC)='OVERV'
      ENDIF
   10 IF(LNEWXS(16)) THEN
        DO IED=1,NEDMAC
          IF(HVECT(IED).EQ.'NFTOT') GO TO 15
        ENDDO
        NEDMAC=NEDMAC+1
        IF(NEDMAC.GT.MAXNED) CALL XABORT('MACPXS: MAXNED OVERFLOW(4).')
        HVECT(NEDMAC)='NFTOT'
      ENDIF
   15 IF(NEDMAC.GT.0) CALL LCMPTC(IPLIST,'ADDXSNAME-P0',8,NEDMAC,HVECT)
*
      JPLIST=LCMLID(IPLIST,'GROUP',NGROUP)
      DO 100 IGR=1,NGROUP
        KPLIST=LCMDIL(JPLIST,IGR)
*----
*  PUT TOTAL, FIXE AND TRANC ON THE MACROLIB IF MODIFIED
*----
        IF(LNEWXS(1)) CALL LCMPUT(KPLIST,'NTOT0',NBMIX,2,XSTOTL(1,IGR))
        IF(LNEWXS(3)) CALL LCMPUT(KPLIST,'FIXE' ,NBMIX,2,XSFIXE(1,IGR))
        IF(LNEWXS(6)) CALL LCMPUT(KPLIST,'TRANC',NBMIX,2,XSTRAN(1,IGR))
        IF(LNEWXS(7)) CALL LCMPUT(KPLIST,'DIFF',NBMIX,2,XSDIFF(1,IGR))
        IF(LNEWXS(8)) CALL LCMPUT(KPLIST,'H-FACTOR',NBMIX,2,XSH(1,IGR))
        IF(LNEWXS(9)) CALL LCMPUT(KPLIST,'NTOT1',NBMIX,2,XSTOT1(1,IGR))
*----
*  PUT CHI AND NUSIGF ON THE MACROLIB IF MODIFIED
*----
        IF(LNEWXS(2))
     >   CALL LCMPUT(KPLIST,'NUSIGF',NBMIX*NIFISS,2,XSFISS(1,1,IGR))
        IF(LNEWXS(4))
     >   CALL LCMPUT(KPLIST,'CHI',NBMIX*NIFISS,2,XSSPEC(1,1,IGR))
*----
*  PUT DIFFX, DIFFY AND DIFFZ ON THE MACROLIB IF MODIFIED
*----
        IF(LNEWXS(10))CALL LCMPUT(KPLIST,'DIFFX',NBMIX,2,XSDIFX(1,IGR))
        IF(LNEWXS(11))CALL LCMPUT(KPLIST,'DIFFY',NBMIX,2,XSDIFY(1,IGR))
        IF(LNEWXS(12))CALL LCMPUT(KPLIST,'DIFFZ',NBMIX,2,XSDIFZ(1,IGR))
*----
*  PUT CHID, NUSIGD AND OVERV ON THE MACROLIB IF MODIFIED
*----
        IF(LNEWXS(13)) THEN
          DO I=1,NDELG
            WRITE(NUSIGD,'(A6,I2.2)') 'NUSIGF',I
            CALL LCMPUT(KPLIST,NUSIGD,NBMIX*NIFISS,2,XSNUDL(1,1,I,IGR))
          ENDDO
        ENDIF
        IF(LNEWXS(14)) THEN
          DO I=1,NDELG
            WRITE(CHID,'(A3,I2.2)') 'CHI',I
            CALL LCMPUT(KPLIST,CHID,NBMIX*NIFISS,2,XSCHDL(1,1,I,IGR))
          ENDDO
        ENDIF
        IF(LNEWXS(15))CALL LCMPUT(KPLIST,'OVERV',NBMIX,2,XSOVRV(1,IGR))
        IF(LNEWXS(16))CALL LCMPUT(KPLIST,'NFTOT',NBMIX,2,XSNFTO(1,IGR))
        IF(LNEWXS(17))CALL LCMPUT(KPLIST,'FLUX-INTG',NBMIX,2,
     >                XSINT0(1,IGR))
        IF(LNEWXS(18))CALL LCMPUT(KPLIST,'FLUX-INTG-P1',NBMIX,2,
     >                XSINT1(1,IGR))
*----
*  COMPRESS AND PUT ON SCATT ON THE MACROLIB IF MODIFIED
*----
        DO 60 IANIS=1,NANISO
          WRITE(CANISO,'(I2.2)') IANIS-1
          IF(LNEWXS(5).AND.ISCATA(IANIS).EQ.2) THEN
            NELEM=0
            DO 50 INM=1,NBMIX
              J2=IGR
              J1=IGR
              DO 20 JGR=1,NGROUP
              IF(XSSCAT(JGR,INM,IANIS,IGR).NE.0.0) THEN
                J2=MAX(J2,JGR)
                J1=MIN(J1,JGR)
              ENDIF
  20          CONTINUE
              INGSCT(INM)=J2-J1+1
              IFGSCT(INM)=J2
              IPOSCT(INM)=NELEM+1
              DO 30 JGR=J2,J1,-1
                NELEM=NELEM+1
                XSWORK(2*NBMIX+NELEM)=XSSCAT(JGR,INM,IANIS,IGR)
  30          CONTINUE
*----
*  STORE DIAGONAL ELEMENTS OF SCATTERING MATRIX
*  AND TOTAL SCATTERING OUT OF GROUP
*----
              XSWORK(INM)=XSSCAT(IGR,INM,IANIS,IGR)
              XSTOT=0.0
              DO 40 JGR=1,NGROUP
              XSTOT=XSTOT+XSSCAT(IGR,INM,IANIS,JGR)
  40          CONTINUE
              XSWORK(NBMIX+INM)=XSTOT
  50        CONTINUE
            CALL LCMPUT(KPLIST,'NJJS'//CANISO,NBMIX,1,INGSCT)
            CALL LCMPUT(KPLIST,'IJJS'//CANISO,NBMIX,1,IFGSCT)
            CALL LCMPUT(KPLIST,'IPOS'//CANISO,NBMIX,1,IPOSCT)
            CALL LCMPUT(KPLIST,'SIGW'//CANISO,NBMIX,2,XSWORK)
            CALL LCMPUT(KPLIST,'SIGS'//CANISO,NBMIX,2,XSWORK(NBMIX+1))
            CALL LCMPUT(KPLIST,'SCAT'//CANISO,NELEM,2,
     >      XSWORK(2*NBMIX+1))
          ENDIF
  60    CONTINUE
*----
*  COMPUTE/RECOVER TRANSPORT CORRECTION
*----
        IF(ITRANC.EQ.2) THEN
*         RECOVER TRANSPORT CORRECTION FROM RECORD 'TRANC'.
          CALL LCMLEN(KPLIST,'TRANC',ILCMLN,ITYLCM)
          IF(ILCMLN.NE.NBMIX) CALL XABORT('MACPXS: NO TRANC RECORD '
     >    //'AVAILABLE')
        ELSE IF(ITRANC.NE.0) THEN
          CALL XDRSET(XSWORK,NBMIX,0.0)
          CALL LCMLEN(KPLIST,'NTOT1',ILCMLN,ITYLCM)
          IF(ILCMLN.EQ.NBMIX) THEN
*           LEAKAGE CORRECTION.
            CALL LCMGET(KPLIST,'NTOT0',XSWORK(NBMIX+1))
            CALL LCMGET(KPLIST,'NTOT1',XSWORK(2*NBMIX+1))
            DO 70 INM=1,NBMIX
            XSWORK(INM)=XSWORK(NBMIX+INM)-XSWORK(2*NBMIX+INM)
  70        CONTINUE
          ENDIF
          IF(ITRANC.EQ.1) THEN
*           APOLLO-TYPE TRANSPORT CORRECTION.
            CALL LCMLEN(KPLIST,'SIGS01',ILCMLN,ITYLCM)
            IF(ILCMLN.EQ.NBMIX) THEN
              CALL LCMGET(KPLIST,'SIGS01',XSWORK(NBMIX+1))
              DO 80 INM=1,NBMIX
              XSWORK(INM)=XSWORK(INM)+XSWORK(NBMIX+INM)
  80          CONTINUE
            ENDIF
          ELSE IF(ITRANC.NE.4) THEN
            CALL XABORT('MACPXS: UNKNOWN TYPE OF CORRECTION.')
          ENDIF
*         ***CAUTION*** 'TRANC' CONTAINS BOTH TRANSPORT AND LEAKAGE
*         CORRECTIONS.
          CALL LCMPUT(KPLIST,'TRANC',NBMIX,2,XSWORK)
        ENDIF
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSWORK)
      DEALLOCATE(IPOSCT,IFGSCT,INGSCT)
      RETURN
      END
