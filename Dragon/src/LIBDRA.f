*DECK LIBDRA
      SUBROUTINE LIBDRA (IPLIB,IPDRL,NAMFIL,NGRO,NBISO,NL,ISONAM,
     1 ISONRF,IPISO,TN,SN,SB,MASKI,NED,HVECT,IMPX,NGF,NGFR,NDEL,NBESP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from a microscopic x-section library (draglib format) to LCM
* data structures.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPDRL   pointer to the draglib (L_DRAGLIB signature).
* NAMFIL  name of the Dragon library file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* TN      temperature of each isotope.
* SN      dilution cross section in each energy group of each
*         isotope. a value of 1.0e10 is used for infinite dilution.
* SB      dilution cross section as used in livolant and jeanpierre
*         normalization.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NDEL    number of precursor groups for delayed neutrons.
* NBESP   number of energy-dependent fission spectra.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXDEL=10,MAXESP=4)
      CHARACTER*(*) HVECT(NED),NAMFIL
      TYPE(C_PTR) IPLIB,IPDRL,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),NED,IMPX,
     1 NGF,NGFR,NDEL,NBESP
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER CD*4,HSMG*131,HTITLE*80,HNISOR*12,HNAMIS*12,HNUSIG*12,
     1 HCHI*12
      PARAMETER (IOUT=6,MAXTMP=50,NOTX=3)
      TYPE(C_PTR) KPLIB
      LOGICAL LSIGF,LGOLD,LOGT,LNZERO
      INTEGER IESP(MAXESP+1)
      DOUBLE PRECISION FACTOR,TERP(MAXTMP)
      REAL TEMP(MAXTMP),ZLAMB(MAXDEL),EESP(MAXESP+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NFS,ITYPRO,ITITLE
      REAL, ALLOCATABLE, DIMENSION(:) :: AWR,DELTA,TOTAL,GOLD,ZNPHI,
     1 ENER,BIN,EBIN,SIGS2,SCAT2,TOTAL2,SIGF2,CHI2,SADD2,GOLD2,BIN2,
     2 ZNPHI2,CHI4G2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS,SIGF,CHI,SADD,CHI4G
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LSCAT,LADD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NFS(NGRO),ITYPRO(NL))
      ALLOCATE(AWR(NBISO),DELTA(NGRO),SIGS(NGRO,NL),SCAT(NGRO,NGRO,NL),
     1 TOTAL(NGRO),SIGF(NGRO,0:MAXDEL),CHI(NGRO,0:MAXDEL),
     2 SADD(NGRO,NED),GOLD(NGRO),ZNPHI(NGRO))
      ALLOCATE(LSCAT(NL),LADD(NED))
*----
*  RECOVER THE GROUP STRUCTURE.
*----
      NGF=NGRO+1
      NGFR=0
      NDEL=0
      IF(IMPX.GT.0) WRITE (IOUT,900) NAMFIL
      CALL LCMLEN(IPDRL,'README',LENGT,ITYLCM)
      IF((IMPX.GT.0).AND.(LENGT.GT.0)) THEN
         ALLOCATE(ITITLE(LENGT))
         CALL LCMGET(IPDRL,'README',ITITLE)
         WRITE (IOUT,940)
         I2=0
         DO 10 J=0,LENGT/20
         I1=I2+1
         I2=MIN(I1+19,LENGT)
         WRITE (HTITLE,'(20A4)') (ITITLE(I),I=I1,I2)
         WRITE (IOUT,'(1X,A80)') HTITLE
   10    CONTINUE
         DEALLOCATE(ITITLE)
      ENDIF
      ALLOCATE(ENER(NGRO+1))
      CALL LCMLEN(IPDRL,'ENERGY',LENGT,ITYLCM)
      LENGT=LENGT-1
      IF(LENGT.NE.NGRO) CALL XABORT('LIBDRA: INVALID GROUP STRUCTURE.')
      CALL LCMGET(IPDRL,'ENERGY',ENER)
      CALL LCMLEN(IPDRL,'DELTAU',LENGT,ITYLCM)
      IF(LENGT.EQ.NGRO) THEN
         CALL LCMGET(IPDRL,'DELTAU',DELTA)
      ELSE IF(LENGT.EQ.0) THEN
         IF(ENER(NGRO+1).EQ.0.0) ENER(NGRO+1)=1.0E-5
         DO 15 J=1,NGRO
         DELTA(J)=LOG(ENER(J)/ENER(J+1))
   15    CONTINUE
      ENDIF
      CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENER)
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
      DEALLOCATE(ENER)
      CALL LCMLEN(IPDRL,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('LIBDRA: MAXESP OVERFLOW.')
         CALL LCMGET(IPDRL,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLIB,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPDRL,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLIB,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
      ALLOCATE(CHI4G(NGRO,NBESP))
*----
*  READ THROUGH DRAGON FILE AND ACCUMULATE CROSS SECTIONS FOR THIS RANGE
*  OF MATS, LEGENDRE ORDERS, AND GROUPS.
*----
      DO 400 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMX),I0=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         CALL LCMLEN(IPDRL,HNISOR,LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL LCMLIB(IPDRL)
            WRITE (HSMG,910) HNAMIS,HNISOR,NAMFIL
            CALL XABORT(HSMG)
         ENDIF
         IF(IMPX.GT.0) WRITE (IOUT,920) HNAMIS,HNISOR
         CALL LCMSIX(IPDRL,HNISOR,1)
*
         CALL LCMGET(IPDRL,'AWR',AWR(IMX))
         CALL LCMLEN(IPDRL,'README',LTITLE,ITYLCM)
         IF(LTITLE.GT.0) THEN
            ALLOCATE(ITITLE(LTITLE))
            CALL LCMGET(IPDRL,'README',ITITLE)
            IF(IMPX.GT.0) THEN
               WRITE (IOUT,930)
               I2=0
               DO 20 J=0,LTITLE/20
               I1=I2+1
               I2=MIN(I1+19,LTITLE)
               WRITE (HTITLE,'(20A4)') (ITITLE(I),I=I1,I2)
               WRITE (IOUT,'(1X,A80)') HTITLE
   20          CONTINUE
            ENDIF
         ENDIF
*----
*  RECOVER BIN TYPE INFORMATION (IF AVAILABLE).
*----
         LBIN=0
         CALL LCMLEN (IPDRL,'BIN-NFS',LENGT,ITYXSM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET (IPDRL,'BIN-NFS',NFS)
            DO 30 I=1,NGRO
            LBIN=LBIN+NFS(I)
   30       CONTINUE
            ALLOCATE(BIN(3*LBIN),EBIN(LBIN+1))
            CALL LCMGET (IPDRL,'BIN-ENERGY',EBIN)
         ENDIF
*
         CALL LCMLEN (IPDRL,'TEMPERATURE',NTMP,ITYLCM)
         IF(NTMP.GT.MAXTMP) CALL XABORT('LIBDRA: MAXTMP OVERFLOW.')
         IF(NTMP.EQ.0) THEN
            CALL LCMLEN (IPDRL,'LAMBDA-D',NDEL0,ITYLCM)
            NDEL=MAX(NDEL,NDEL0)
            IF(NDEL0.GT.MAXDEL) CALL XABORT('LIBDRA: MAXDEL OVERFLOW.')
            IF(NDEL0.GT.0) CALL LCMGET (IPDRL,'LAMBDA-D',ZLAMB)
            CALL LIBDRB (IPDRL,NGRO,NL,NDEL0,NBESP,SN(1,IMX),SB(1,IMX),
     1      NED,HVECT,DELTA,LBIN,NFS,EBIN,IMPX,NGF,NGFR,LSCAT,LSIGF,
     2      LADD,LGOLD,SIGS,SCAT,TOTAL,ZNPHI,SIGF,CHI,CHI4G,SADD,GOLD,
     3      BIN)
         ELSE
*----
*  PERFORM TEMPERATURE LAGRANGIAN INTERPOLATION (ORDER ABS(NOTX)).
*----
            CALL LCMSIX (IPDRL,'SUBTMP0001',1)
            CALL LCMLEN (IPDRL,'LAMBDA-D',NDEL0,ITYLCM)
            NDEL=MAX(NDEL,NDEL0)
            IF(NDEL0.GT.MAXDEL) CALL XABORT('LIBDRA: MAXDEL OVERFLOW.')
            IF(NDEL0.GT.0) CALL LCMGET (IPDRL,'LAMBDA-D',ZLAMB)
            CALL LCMSIX (IPDRL,' ',2)
            CALL LCMGET (IPDRL,'TEMPERATURE',TEMP)
            CALL LIBLEX(NTMP,TN(IMX),TEMP,NOTX,TERP)
            DO 121 IG1=1,NGRO
            TOTAL(IG1)=0.0
            ZNPHI(IG1)=0.0
            DO 100 IDEL=0,NDEL0
            SIGF(IG1,IDEL)=0.0
            CHI(IG1,IDEL)=0.0
  100       CONTINUE
            DO 105 ISP=1,NBESP
            CHI4G(IG1,ISP)=0.0
  105       CONTINUE
            GOLD(IG1)=0.0
            DO 115 IL=1,NL
            SIGS(IG1,IL)=0.0
            DO 110 IG2=1,NGRO
            SCAT(IG1,IG2,IL)=0.0
  110       CONTINUE
  115       CONTINUE
            DO 120 IED=1,NED
            SADD(IG1,IED)=0.0
  120       CONTINUE
  121       CONTINUE
            DO 125 IG=1,3*LBIN
            BIN(IG)=0.0
  125       CONTINUE
            ALLOCATE(SIGS2(NGRO*NL),SCAT2(NGRO*NGRO*NL),TOTAL2(NGRO),
     1      SIGF2(NGRO*(NDEL0+1)),CHI2(NGRO*(NDEL0+1)),SADD2(NGRO*NED),
     2      GOLD2(NGRO),BIN2(3*LBIN),ZNPHI2(NGRO),CHI4G2(NGRO*NBESP))
            FACTOR=1.0D0
            DO 210 ITM=1,NTMP
            TERPM=REAL(TERP(ITM))
            FACTOR=FACTOR-TERP(ITM)
            IF(TERPM.EQ.0.0) GO TO 210
            IF(IMPX.GT.4) WRITE(6,'(/30H DRAGLIB ACCESS AT TEMPERATURE,
     >      1P,E12.4,18H KELVIN. FACTOR = ,E12.4)') TEMP(ITM),TERPM
            WRITE (CD,'(I4.4)') ITM
            CALL LCMSIX (IPDRL,'SUBTMP'//CD,1)
            CALL LIBDRB (IPDRL,NGRO,NL,NDEL0,NBESP,SN(1,IMX),SB(1,IMX),
     1      NED,HVECT,DELTA,LBIN,NFS,EBIN,IMPX,NGF,NGFR,LSCAT,LSIGF,
     2      LADD,LGOLD,SIGS2,SCAT2,TOTAL2,ZNPHI2,SIGF2,CHI2,CHI4G2,
     3      SADD2,GOLD2,BIN2)
            CALL LCMSIX (IPDRL,' ',2)
            DO 130 IG=1,NGRO
              TOTAL(IG)=TOTAL(IG)+TERPM*TOTAL2(IG)
              ZNPHI(IG)=ZNPHI(IG)+TERPM*ZNPHI2(IG)
  130       CONTINUE
            IF(LSIGF) THEN
              DO 141 IDEL=0,NDEL0
              DO 140 IG=1,NGRO
                IOFSET=IDEL*NGRO+IG-1
                SIGF(IG,IDEL)=SIGF(IG,IDEL)+TERPM*SIGF2(IOFSET+1)
                CHI(IG,IDEL)=CHI(IG,IDEL)+TERPM*CHI2(IOFSET+1)
  140         CONTINUE
  141         CONTINUE
              DO 146 ISP=1,NBESP
              DO 145 IG=1,NGRO
                IOFSET=(ISP-1)*NGRO+IG-1
                CHI4G(IG,ISP)=CHI4G(IG,ISP)+TERPM*CHI4G2(IOFSET+1)
  145         CONTINUE
  146         CONTINUE
            ENDIF
            DO 160 IL=1,NL
              IF(LSCAT(IL)) THEN
                DO 150 IG2=1,NGRO
                  SIGS(IG2,IL)=SIGS(IG2,IL)+TERPM*SIGS2((IL-1)*NGRO+IG2)
                  IOF=(IL-1)*NGRO*NGRO+(IG2-1)*NGRO
                  DO 151 IG1=1,NGRO
                    SCAT(IG1,IG2,IL)=SCAT(IG1,IG2,IL)+TERPM*
     >              SCAT2(IOF+IG1)
  151             CONTINUE
  150           CONTINUE
              ENDIF
  160       CONTINUE
            DO 180 IED=1,NED
            IF(LADD(IED)) THEN
               DO 170 IG=1,NGRO
                 SADD(IG,IED)=SADD(IG,IED)+TERPM*SADD2((IED-1)*NGRO+IG)
  170          CONTINUE
            ENDIF
  180       CONTINUE
            IF(LGOLD) THEN
               DO 190 IG=1,NGRO
               GOLD(IG)=GOLD(IG)+TERPM*GOLD2(IG)
  190          CONTINUE
            ENDIF
            DO 200 IG=1,3*LBIN
            BIN(IG)=BIN(IG)+TERPM*BIN2(IG)
  200       CONTINUE
  210       CONTINUE
            DEALLOCATE(CHI4G2,ZNPHI2,BIN2,GOLD2,SADD2,CHI2,SIGF2,TOTAL2,
     >      SCAT2,SIGS2)
            IF(ABS(FACTOR).GT.1.0D-4) CALL XABORT('LIBDRA: TERP ERROR')
         ENDIF
         CALL LCMSIX(IPDRL,' ',2)
*----
*        SAVE CROSS SECTION DATA ON LCM.
*----
         KPLIB=IPISO(IMX) ! set IMX-th isotope
         CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
         CALL LCMPUT(KPLIB,'AWR',1,2,AWR(IMX))
         IF(LTITLE.GT.0) THEN
            CALL LCMPUT(KPLIB,'README',LTITLE,3,ITITLE)
            DEALLOCATE(ITITLE)
         ENDIF
         DO 220 IG=1,NGRO
         IF(TOTAL(IG).LT.0.0) THEN
            WRITE(HSMG,'(42HLIBDRA: NEGATIVE TOTAL CROSS SECTION IN GR,
     1      3HOUP,I4,14H FOR ISOTOPE '',A12,2H''.)') IG,HNAMIS
            CALL XABORT(HSMG)
         ELSE IF(ZNPHI(IG).LT.0.0) THEN
            WRITE(HSMG,'(41HLIBDRA: NEGATIVE INTEGRATED FLUX IN GROUP,
     1      I4,14H FOR ISOTOPE '',A12,2H''.)') IG,HNAMIS
            CALL XABORT(HSMG)
         ENDIF
  220    CONTINUE
         CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,TOTAL)
         CALL LCMPUT(KPLIB,'NWT0',NGRO,2,ZNPHI)
         IF(NDEL0.GT.0) CALL LCMPUT (KPLIB,'LAMBDA-D',NDEL0,2,ZLAMB)
         IF(LSIGF) THEN
            DO 250 IDEL=0,NDEL0
            IF(IDEL.EQ.0) THEN
               HNUSIG='NUSIGF'
            ELSE
               WRITE(HNUSIG,'(6HNUSIGF,I2.2)') IDEL
            ENDIF
            CALL LCMPUT(KPLIB,HNUSIG,NGRO,2,SIGF(1,IDEL))
            IF(IDEL.EQ.0) THEN
               IF(NBESP.GT.0) GO TO 250
               HCHI='CHI'
            ELSE
               WRITE(HCHI,'(3HCHI,I2.2)') IDEL
            ENDIF
            CALL LCMPUT(KPLIB,HCHI,NGRO,2,CHI(1,IDEL))
  250       CONTINUE
            DO 260 ISP=1,NBESP
               LNZERO=.FALSE.
               DO 255 IG=1,NGRO
               LNZERO=LNZERO.OR.(CHI4G(IG,ISP).NE.0.0)
  255          CONTINUE
               IF(LNZERO) THEN
                 WRITE(HCHI,'(5HCHI--,I2.2)') ISP
                 CALL LCMPUT(KPLIB,HCHI,NGRO,2,CHI4G(1,ISP))
               ENDIF
  260       CONTINUE
         ENDIF
         CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
         DO 340 IED=1,NED
         IF(LADD(IED).AND.(HVECT(IED)(:3).NE.'CHI')
     1               .AND.(HVECT(IED)(:2).NE.'NU')
     2               .AND.(HVECT(IED).NE.'NTOT0')
     3               .AND.(HVECT(IED)(:3).NE.'NWT')) THEN
            CALL LCMPUT(KPLIB,HVECT(IED),NGRO,2,SADD(1,IED))
         ENDIF
  340    CONTINUE
         IF(LGOLD) CALL LCMPUT(KPLIB,'NGOLD',NGRO,2,GOLD)
         IF(LBIN.GT.0) THEN
            CALL LCMPUT(KPLIB,'BIN-NFS',NGRO,1,NFS)
            CALL LCMPUT(KPLIB,'BIN-ENERGY',LBIN+1,2,EBIN)
            CALL LCMPUT(KPLIB,'BIN-NTOT0',LBIN,2,BIN)
            CALL LCMPUT(KPLIB,'BIN-SIGS00',LBIN,2,BIN(LBIN+1))
            LOGT=.FALSE.
            DO 350 I=1,LBIN
            LOGT=LOGT.OR.(BIN(2*LBIN+I).NE.0.0)
  350       CONTINUE
            IF(LOGT) THEN
               CALL LCMPUT(KPLIB,'BIN-NUSIGF',LBIN,2,BIN(2*LBIN+1))
            ENDIF
            DEALLOCATE(EBIN,BIN)
         ENDIF
         IF(IMPX.GT.9) CALL LCMLIB(KPLIB)
      ENDIF
  400 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LADD,LSCAT)
      DEALLOCATE(CHI4G,ZNPHI,GOLD,SADD,CHI,SIGF,TOTAL,SCAT,SIGS,DELTA,
     1 AWR)
      DEALLOCATE(ITYPRO,NFS)
      RETURN
*
  900 FORMAT(/33H PROCESSING DRAGON LIBRARY NAMED ,A12,1H.)
  910 FORMAT(26HLIBDRA: MATERIAL/ISOTOPE ',A12,5H' = ',A12,9H' IS MISS,
     1 25HING ON DRAGON FILE NAMED ,A12,1H.)
  920 FORMAT(/30H PROCESSING ISOTOPE/MATERIAL ',A12,11H' (HNISOR=',A12,
     1 3H').)
  930 FORMAT(/23H ISOTOPE/MATERIAL INFO:)
  940 FORMAT(/24H X-SECTION LIBRARY INFO:)
      END
