*DECK MACWXS
      SUBROUTINE MACWXS(IPLIST,IPRINT,NGROUP,NBMIX,NIFISS,NANISO,ICTRA,
     > NEDMAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print/check cross section information if required.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIST  LCM pointer to the macrolib.
* IPRINT  print level.
* NGROUP  number of energy groups.
* NBMIX   number of mixtures.
* NIFISS  number of fissile isotopes.
* NANISO  maximum Legendre order:
*         =1 isotropic collision;
*         =2 linearly anisotropic collision.
* ICTRA   type of transport correction.
* NEDMAC  number of macro edit cross sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER    IPRINT,NGROUP,NBMIX,NIFISS,NANISO,ICTRA,NEDMAC
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIST,KPLIST
      CHARACTER  CANISO*2,NAMREC*12,CEDNAM*12,NAMRE1*12,NAMRE2*12
      PARAMETER (IUNOUT=6,ILCMUP=1,ILCMDN=2,IOUT=6)
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGSCT,IFGSCT,IPOSCT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CADNAM
      REAL, ALLOCATABLE, DIMENSION(:) :: XSWORK,XSWOR2,TCOR
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCWORK
*----
*  SCRATCH STORAGE ALLOCATION
*   INGSCT  number of scattering group for cross sections.
*   IFGSCT  first scattering group for cross sections.
*   IPOSCT  mixture location in cross section matrix.
*   TCOR    transport correction.
*   CADNAM  additional cross section names.
*----
      ALLOCATE(INGSCT(NBMIX),IFGSCT(NBMIX),IPOSCT(NBMIX),
     > CADNAM(2,NEDMAC+1))
      ALLOCATE(XSWORK(NBMIX*(2+NGROUP)),XSWOR2(NBMIX*NIFISS),
     > SCWORK(NBMIX,NANISO,NGROUP),TCOR(NBMIX))
*
      IF(IPRINT.GT.1) WRITE(IUNOUT,6000)
      CALL XDRSET(SCWORK,NBMIX*NANISO*NGROUP,0.0)
      IF(NEDMAC.GT.0) CALL LCMGET(IPLIST,'ADDXSNAME-P0',CADNAM)
      JPLIST=LCMGID(IPLIST,'GROUP')
      DO 100 IGR=1,NGROUP
        KPLIST=LCMGIL(JPLIST,IGR)
        IF(ICTRA.EQ.0) THEN
          CALL XDRSET(TCOR,NBMIX,0.0)
        ELSE
          CALL LCMGET(KPLIST,'TRANC',TCOR)
        ENDIF
*----
*  PRINT TOTAL X-S AND FIXE SOURCES INFORMATION
*----
        IF(IPRINT.GT.1) THEN
          IF(NBMIX.LE.7) THEN
            WRITE(IUNOUT,6001) IGR,(II,II=1,NBMIX)
          ELSE
            WRITE(IUNOUT,6001) IGR,(II,II=1,7)
            WRITE(IUNOUT,6011) (II,II=8,NBMIX)
          ENDIF
          CALL LCMGET(KPLIST,'NTOT0',XSWORK)
          IF(NBMIX.LE.7) THEN
            WRITE(IUNOUT,6002) 'NTOT0       ',
     >      (XSWORK(II)-TCOR(II),II=1,NBMIX)
          ELSE
            WRITE(IUNOUT,6002) 'NTOT0       ',
     >                         (XSWORK(II)-TCOR(II),II=1,7)
            WRITE(IUNOUT,6010) (XSWORK(II)-TCOR(II),II=8,NBMIX)
          ENDIF
          CALL LCMLEN(KPLIST,'NTOT1',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'NTOT1',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'NTOT1       ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'NTOT1       ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'FIXE',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'FIXE',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'FIXE        ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'FIXE        ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'TRANC',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'TRANC',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'TRANC       ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'TRANC       ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'DIFF',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'DIFF',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'DIFF        ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'DIFF        ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'H-FACTOR',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'H-FACTOR',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'H-FACTOR    ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'H-FACTOR    ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'FLUX-INTG',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'FLUX-INTG',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'FLUX-INTG   ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'FLUX-INTG   ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
          CALL LCMLEN(KPLIST,'FLUX-INTG-P1',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'FLUX-INTG-P1',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'FLUX-INTG-P1',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'FLUX-INTG-P1',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
        ENDIF
*----
*  PRINT FISSION INFORMATION
*----
        IF(IPRINT.GT.1) THEN
          CALL LCMLEN(KPLIST,'NUSIGF',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'NUSIGF',XSWOR2)
            DO 110 IFISS=1,NIFISS
              IOF=(IFISS-1)*NBMIX
              WRITE(IUNOUT,6003) IFISS
              IF(NBMIX.LE.7) THEN
                WRITE(IUNOUT,6002) 'NUSIGF',(XSWOR2(IOF+II),II=1,NBMIX)
              ELSE
                WRITE(IUNOUT,6002) 'NUSIGF',(XSWOR2(IOF+II),II=1,7)
                WRITE(IUNOUT,6010) (XSWOR2(IOF+II),II=8,NBMIX)
              ENDIF
 110        CONTINUE
          ENDIF
*
          CALL LCMLEN(KPLIST,'CHI',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'CHI',XSWOR2)
            DO 115 IFISS=1,NIFISS
              IOF=(IFISS-1)*NBMIX
              WRITE(IUNOUT,6003) IFISS
              IF(NBMIX.LE.7) THEN
                WRITE(IUNOUT,6002) 'CHI',(XSWOR2(IOF+II),II=1,NBMIX)
              ELSE
                WRITE(IUNOUT,6002) 'CHI',(XSWOR2(IOF+II),II=1,7)
                WRITE(IUNOUT,6010) (XSWOR2(IOF+II),II=8,NBMIX)
              ENDIF
 115        CONTINUE
          ENDIF
*
          CALL LCMLEN(KPLIST,'NFTOT',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPLIST,'NFTOT',XSWORK)
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) 'NFTOT    ',(XSWORK(II),II=1,NBMIX)
            ELSE
              WRITE(IUNOUT,6002) 'NFTOT    ',(XSWORK(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II),II=8,NBMIX)
            ENDIF
          ENDIF
        ENDIF
*----
*  PRINT SCATTERING MATRICES INFORMATION
*----
        MAXMXR=0
        IF(ICTRA.NE.0) THEN
          NNANIS=1
        ELSE
          NNANIS=NANISO
        ENDIF
        DO 120 IANIS=1,NNANIS
          WRITE(CANISO,'(I2.2)') IANIS-1
          CALL LCMLEN(KPLIST,'NJJS'//CANISO,ILCMLN,ITYLCM)
          IF(ILCMLN.EQ.0) GO TO 120
          IF(ILCMLN.LT.NBMIX) THEN
            MAXMXR=ILCMLN
          ELSE IF(ILCMLN.GT.NBMIX) THEN
            CALL XABORT('MACWXS: NUMBER OF MIXTURES ON LCM IS LARGER'//
     >                  ' THAN THE MAXIMUM PROVIDED ON NMIX CARD')
          ELSE
            MAXMXR=NBMIX
          ENDIF
          NAMREC='SIGW'//CANISO
          NAMRE1='SIGS'//CANISO
          NAMRE2='SCAT'//CANISO
          CALL LCMGET(KPLIST,'NJJS'//CANISO,INGSCT)
          CALL LCMGET(KPLIST,'IJJS'//CANISO,IFGSCT)
          CALL LCMGET(KPLIST,'IPOS'//CANISO,IPOSCT)
          CALL LCMGET(KPLIST,NAMREC,XSWORK)
          CALL LCMLEN(KPLIST,NAMRE1,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) CALL LCMGET(KPLIST,NAMRE1,XSWORK(NBMIX+1))
          CALL LCMGET(KPLIST,NAMRE2,XSWORK(2*NBMIX+1))
          IF(IPRINT.GT.1) THEN
            IF(NBMIX.LE.7) THEN
              WRITE(IUNOUT,6002) NAMREC,
     >          (XSWORK(II)-TCOR(II),II=1,MAXMXR)
              IF(ILCMLN.GT.0) THEN
                WRITE(IUNOUT,6002) NAMRE1,
     >            (XSWORK(NBMIX+II)-TCOR(II),II=1,MAXMXR)
              ENDIF
            ELSE
              WRITE(IUNOUT,6002) NAMREC,
     >          (XSWORK(II)-TCOR(II),II=1,7)
              WRITE(IUNOUT,6010) (XSWORK(II)-TCOR(II),II=8,MAXMXR)
              IF(ILCMLN.GT.0) THEN
                WRITE(IUNOUT,6002) NAMRE1,
     >            (XSWORK(NBMIX+II)-TCOR(II),II=1,7)
                WRITE(IUNOUT,6010) (XSWORK(NBMIX+II)-TCOR(II),II=8,
     >            MAXMXR)
              ENDIF
            ENDIF
          ENDIF
          IF(IPRINT.GT.2) THEN
            IF=2*NBMIX+1
            IS=2*NBMIX+INGSCT(1)
            IDG=IFGSCT(1)
            IFG=IDG-INGSCT(1)+1
            IDIAG=IS+IFG-IGR
            XSWORK(IDIAG)=XSWORK(IDIAG)-TCOR(1)
            IF(INGSCT(1).GT.0) THEN
              IF(IS-IF+1.LE.4) THEN
                WRITE(IUNOUT,6004) NAMRE2,1,IFG,IDG,
     >                            (XSWORK(J),J=IS,IF,-1)
              ELSE
                WRITE(IUNOUT,6004) NAMRE2,1,IFG,IDG,
     >                            (XSWORK(J),J=IS,IS-3,-1)
                WRITE(IUNOUT,6010)(XSWORK(J),J=IS-4,IF,-1)
              ENDIF
            ENDIF
            DO 130 IMAT=2,MAXMXR
              IF=IS+1
              IS=IS+INGSCT(IMAT)
              IDG=IFGSCT(IMAT)
              IFG=IDG-INGSCT(IMAT)+1
              IDIAG=IS+IFG-IGR
              XSWORK(IDIAG)=XSWORK(IDIAG)-TCOR(IMAT)
              IF(INGSCT(IMAT).GT.0) THEN
                IF(IS-IF+1.LE.4) THEN
                  WRITE(IUNOUT,6005) IMAT,IFG,IDG,
     >                              (XSWORK(J),J=IS,IF,-1)
                ELSE
                  WRITE(IUNOUT,6005) IMAT,IFG,IDG,
     >                              (XSWORK(J),J=IS,IS-3,-1)
                  WRITE(IUNOUT,6010)(XSWORK(J),J=IS-4,IF,-1)
                ENDIF
              ENDIF
 130        CONTINUE
          ENDIF
          IF((IPRINT.GT.4).OR.(IPRINT.LT.0)) THEN
*----
*  VALIDATION OF SCATTERING INFORMATION
*----
            IPO=0
            DO 90 IMAT=1,NBMIX
              IDG=IFGSCT(IMAT)
              IFG=IDG-INGSCT(IMAT)+1
              IF(IPOSCT(IMAT).NE.IPO+1) THEN
                WRITE(IUNOUT,6006) 'IPOS'//CANISO,IGR,IMAT
              ENDIF
              XS=XSWORK(IMAT)
              IF(XS.GT.0.0) THEN
                ERR=ABS((XS-XSWORK(2*NBMIX+IPO+1+IDG-IGR))/XS)
                IF(ERR.GT.1.0E-5) THEN
                  WRITE(IUNOUT,6006) 'SIGW'//CANISO,IGR,IMAT
                  WRITE(IUNOUT,'(7H ERROR=,1P,E15.5)') ERR
                ENDIF
              ENDIF
              DO 80 JGR=IDG,IFG,-1
                IPO=IPO+1
                SCWORK(IMAT,IANIS,JGR)=SCWORK(IMAT,IANIS,JGR)
     >          +XSWORK(2*NBMIX+IPO)
  80          CONTINUE
  90        CONTINUE
          ENDIF
 120    CONTINUE
        IF(IPRINT.GE.4) THEN
*----
*  PRINT ADDITIONAL XS INFORMATION IF REQUIRED
*----
          DO 200 IED=1,NEDMAC
            WRITE(CEDNAM,'(2A4)') CADNAM(1,IED),CADNAM(2,IED)
            CALL LCMLEN(KPLIST,CEDNAM,ILCMLN,ILCMTY)
            IF(ILCMLN .GT. 0) THEN
              CALL LCMGET(KPLIST,CEDNAM,XSWORK)
              IF(NBMIX .LE. 6) THEN
                WRITE(IOUT,6007) CEDNAM,(XSWORK(ITC),ITC=1,NBMIX)
              ELSE
                WRITE(IOUT,6007) CEDNAM,(XSWORK(ITC),ITC=1,6)
                WRITE(IOUT,6010) (XSWORK(ITC),ITC=7,NBMIX)
              ENDIF
            ENDIF
 200      CONTINUE
        ENDIF
 100  CONTINUE
      IF((IPRINT.GT.4).OR.(IPRINT.LT.0)) THEN
        DO 160 IGR=1,NGROUP
        KPLIST=LCMGIL(JPLIST,IGR)
          DO 150 IANIS=1,NANISO
            WRITE(CANISO,'(I2.2)') IANIS-1
            CALL LCMLEN(KPLIST,'SIGS'//CANISO,ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPLIST,'SIGS'//CANISO,XSWORK)
              DO 140 IMAT=1,NBMIX
                XS=XSWORK(IMAT)
                IF(XS.GT.0.0) THEN
                  ERR=ABS((SCWORK(IMAT,IANIS,IGR)-XS)/XS)
                  IF(ERR.GT.1.0E-4) THEN
                    WRITE(IUNOUT,6006) 'SIGS'//CANISO,IGR,IMAT
                    WRITE(IUNOUT,'(7H ERROR=,1P,E15.5)') ERR
                  ENDIF
                ENDIF
 140          CONTINUE
            ENDIF
 150      CONTINUE
 160    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TCOR,SCWORK,XSWOR2,XSWORK)
      DEALLOCATE(CADNAM,IPOSCT,IFGSCT,INGSCT)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(//' EDITION OF X-S STORED ON LCM.'/1X,29(1H-))
 6001 FORMAT(/' G R O U P : ',I5/' MIXTURE',I8,6I16)
 6002 FORMAT(1X,A12,1P,7E16.8)
 6003 FORMAT(/14X,'FISSIONABLE ISOTOPE POSITION:',I5)
 6004 FORMAT(1X,A12,1X,'MIXTURE:',I5,3X,'FROM GROUPS:',I5,' TO ',I5,5X,
     >1P,4E16.8)
 6005 FORMAT(14X,'MIXTURE:',I5,3X,'FROM GROUPS:',I5,' TO ',I5,5X,
     >1P,4E16.8)
 6006 FORMAT(/39H MACWXS: INCONSISTENT VALUE OF RECORD ',A12,7H' IN GR,
     >3HOUP,I4,12H AND MIXTURE,I6,1H.)
 6007 FORMAT(1X,'ADDITIONAL XS : ',A12,1P,6E16.8)
 6010 FORMAT(1P,(13X,7E16.8))
 6011 FORMAT(7I16)
      END
