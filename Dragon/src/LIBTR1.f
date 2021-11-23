*DECK LIBTR1
      SUBROUTINE LIBTR1 (IPLIB,NAMFIL,NGRO,NBISO,NL,ISONAM,ISONRF,
     1 IPISO,ICOHNA,IINCNA,NTFG,TN,SN,SB,MASKI,NED,HVECT,ITIME,IMPX,
     2 NGF,NGFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from matxs to LCM data structures. Use matxs format from NJOY-II
* or NJOY89.
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
* NAMFIL  name of the MATXS library file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library reference name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* ICOHNA  hcoh name.
* IINCNA  hinc name.
* NTFG    number of thermal groups where the thermal inelastic
*         correction is applied.
* TN      temperature of each isotope.
* SN      dilution cross section in each energy group of each
*         isotope. A value of 1.0E10 is used for infinite dilution.
* SB      dilution cross section as used by Livolant and Jeanpierre
*         normalization.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* NED     number of extra vector edits from matxs.
* HVECT   matxs names of the extra vector edits.
*          MATXS reserved names:
*          NWT0/NWT1    p0/p1 library weight function;
*          NTOT0/NTOT1  p0/p1 neutron total cross sections;
*          NELAS  neutron elastic scattering cross section;
*          NINEL  neutron inelastic scattering cross section;
*          NG     radiative capture cross section;
*          NFTOT  total fission cross section;
*          NUDEL  number of delayed secondary neutrons (nu-d);
*          NFSLO  nu * slow fission cross section;
*          CHIS/CHID  slow/delayed fission spectrum;
*          NF/NNF/N2NF/N3NF  nu * partial fission cross sections;
*          N2N/N3N/N4N  (n,2n),(n,3n),(n,4n) cross sections.
* ITIME   MATXS type of fission spectrum:
*         =1 steady-state; =2 prompt.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
*
*Reference:
* R. E. Macfarlane, TRANSX-CTR: A code for interfacing matxs cross-
* section libraries to nuclear transport codes for fusion systems
* analysis, Los Alamos National Laboratory, Report LA-9863-MS,
* New Mexico, February 1984.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT CHARACTER*6 (H)
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),
     1 ICOHNA(2,NBISO),IINCNA(2,NBISO),NTFG(NBISO),NED,ITIME,IMPX,
     2 NGF,NGFR
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      LOGICAL MASKI(NBISO)
      CHARACTER NAMFIL*(*),HVECT(NED)*(*)
*----
*  LOCAL VARIABLES
*----
      CHARACTER FORM*4,HSMG*131,HNISOR*12,HINC*6,HCOH*6,README*88,
     1 HNAMIS*12
      PARAMETER (MULT=2,IOUT=6,FORM='(A6)',MAXA=1000)
      TYPE(C_PTR) KPLIB
      LOGICAL LSUBM1,LTHERM,LTIME,LTERP
      DOUBLE PRECISION HA(MAXA/2)
      REAL A(MAXA)
      INTEGER IA(MAXA),IHGAR(22)
      CHARACTER*6 HGAR(18)
      EQUIVALENCE (A(1),IA(1),HA(1))
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPR
      REAL, ALLOCATABLE, DIMENSION(:) :: AWR,CNORM,SNORM,DNORM,SFIS,
     1 SAVE,VECT,GAR,XS,TERP,TEMP,SIGZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: CHI,SIGF,TOTAL,FLUX
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGS,SCAT
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOGIED
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPR(2,NBISO),ITYPRO(NL))
      ALLOCATE(AWR(NBISO),CNORM(NBISO),SNORM(NBISO),DNORM(NBISO),
     1 SFIS(NGRO),SAVE(NGRO),CHI(NGRO,NBISO),SIGS(NGRO,NL,NBISO),
     2 SIGF(NGRO,NBISO),TOTAL(NGRO,NBISO),SCAT(NGRO,NGRO,NL),
     3 FLUX(NGRO,NBISO),VECT(NGRO),GAR(NGRO))
      ALLOCATE(LOGIED(NED,NBISO))
*
      NGF=NGRO+1
      NGFR=0
      DO 20 I=1,NBISO
      IPR(1,I)=0
      IPR(2,I)=0
   20 CONTINUE
      IF(IMPX.GT.0) WRITE (IOUT,890) NAMFIL
      NIN=KDROPN(NAMFIL,2,2,0)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(36HLIBTR1: UNABLE TO OPEN LIBRARY FILE ,A8,1H.)')
     1   NAMFIL
         CALL XABORT(HSMG)
      ENDIF
*----
*  INITIALIZE MATXS LIBRARY
*----
      NWDS=1+3*MULT
      IREC=1
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      WRITE(HN,FORM) HA(1)
      WRITE(HU,FORM) HA(2)
      WRITE(HS,FORM) HA(3)
      IVER=IA(1+3*MULT)
      IF(IMPX.GT.0) WRITE (IOUT,935) HN,HU,HS,IVER
*----
*  FILE CONTROL
*----
      NWDS=3
      IREC=2
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      NHOLL=IA(3)
*----
*  SET HOLLERITH IDENTIFICATION
*----
      NWDS=NHOLL*MULT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(1).')
      IREC=3
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      WRITE(README(9:),'(6H FROM ,12A6)') (HA(I),I=1,MIN(NHOLL,12))
      IF(IMPX.GT.0) WRITE (IOUT,'(1X,12A6)') (HA(I),I=1,MIN(NHOLL,12))
*----
*  FILE DATA
*----
      NWDS=(NPART+NTYPE)*MULT+6*NTYPE+NPART
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(2).')
      IREC=4
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NWC=NPART+NTYPE
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      L2H=(L2-1)/MULT+1
*----
*  CHECK GROUP STRUCTURES
*----
      NEX1=(NPART+NTYPE)*MULT+6*NTYPE
      DO 170 I=1,NPART
      WRITE(HPART,FORM) HA(I)
      NG=IA(NEX1+I)
      IF(((HPART.EQ.'NEUT').OR.(HPART.EQ.'N')).AND.(NG.NE.NGRO))
     1 CALL XABORT('LIBTR1: INCONSISTENT GROUP STRUCTURES.')
      NWDS=IA(NEX1+I)+1
      ALLOCATE(XS(NWDS))
      IREC=IREC+1
*     ------------------------------
      CALL XDREED (NIN,IREC,XS,NWDS)
*     ------------------------------
      IF((HPART.EQ.'NEUT').OR.(HPART.EQ.'N')) THEN
*        ENERGY BOUND IN EACH GROUP (IN EV):
         CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,XS)
         DO 169 J=1,NGRO
         VECT(J)=LOG(XS(J)/XS(J+1))
  169    CONTINUE
         CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,VECT)
      ENDIF
      DEALLOCATE(XS)
  170 CONTINUE
      IRZT=5+NPART
*----
*  READ THROUGH MATXS FILE AND ACCUMULATE CROSS SECTIONS
*  FOR THIS RANGE OF MATS, LEGENDRE ORDERS, AND GROUPS.
*----
      DO 212 KM=1,NBISO
      DO 205 IED=1,NED
      LOGIED(IED,KM)=.FALSE.
  205 CONTINUE
      CNORM(KM)=0.0
      DO 211 KG=1,NGRO
      CHI(KG,KM)=0.0
      SIGF(KG,KM)=0.0
      TOTAL(KG,KM)=0.0
      DO 210 IL=1,NL
      SIGS(KG,IL,KM)=0.0
  210 CONTINUE
  211 CONTINUE
  212 CONTINUE
*----
*  ***DATA TYPE LOOP***
*----
      DO 680 IT=1,NTYPE
      WRITE(HTYPE,FORM) HA(NPART+IT)
      IF(HTYPE.EQ.'NSCAT') THEN
         ITYPE=1
      ELSE IF(HTYPE.EQ.'NTHERM') THEN
         ITYPE=2
      ELSE
         GO TO 680
      ENDIF
      NDEX=(NPART+NTYPE)*MULT+IT
      NMAT=IA(NDEX)
      NDEX=NDEX+NTYPE
      NINP=IA(NDEX)
      NDEX=NDEX+NTYPE
      NING=IA(NDEX)
      NDEX=NDEX+NTYPE
      NOUTP=IA(NDEX)
      NDEX=NDEX+NTYPE
      NOUTG=IA(NDEX)
      NDEX=NDEX+NTYPE
      LOCT=IA(NDEX)
*----
*  DATA TYPE CONTROL
*----
      NWDS=(2+MULT)*NMAT+NINP+NOUTP+1
      IF(L2+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(3).')
      IREC=LOCT+IRZT
*     ---------------------------------
      CALL XDREED (NIN,IREC,A(L2),NWDS)
*     ---------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      LMC=L2+NWDS
      LMCH=L2H+NWDS/MULT
      NSBLK=IA(L2+NMAT*(MULT+2)+NINP+NOUTP)
      IRZM=IREC+1
*----
*  ***MATERIAL/ISOTOPE LOOP***
*----
      DO 670 IM=1,NMAT
      WRITE (HMAT,FORM) HA(L2H-1+IM)
  300 DO 305 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         IMT=IMX
         WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,IMX),ITC=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(ITC,IMX),ITC=1,3)
         WRITE(HCOH,'(A4,A2)') (ICOHNA(ITC,IMX),ITC=1,2)
         WRITE(HINC,'(A4,A2)') (IINCNA(ITC,IMX),ITC=1,2)
         IF(NTFG(IMX).EQ.0) IPR(2,IMX)=1
         IF((HMAT.EQ.HNISOR(:6)).AND.(IPR(ITYPE,IMX).EQ.0)) GO TO 306
      ENDIF
  305 CONTINUE
      GO TO 670
*----
*  MATERIAL CONTROL
*----
  306 IPR(ITYPE,IMT)=1
      KPLIB=IPISO(IMT) ! set IMT-th isotope
      IF(ITYPE.EQ.1) THEN
         DO 227 IL=0,NL-1
         DO 226 IG2=1,NGRO
         DO 225 IG1=1,NGRO
         SCAT(IG1,IG2,IL+1)=0.0
  225    CONTINUE
  226    CONTINUE
  227    CONTINUE
      ELSE
         CALL XDRLGS(KPLIB,-1,0,0,NL-1,1,NGRO,SIGS(1,1,IMT),SCAT,
     1   ITYPRO)
      ENDIF
*
      LOC=L2-1+MULT*NMAT+IM
      NSUBM=IA(LOC)
      LOCM=IA(LOC+NMAT)
      IREC=LOCM+IRZM
      NWDS=MULT+1+6*NSUBM
      IF(LMC+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(4).')
*     ----------------------------------
      CALL XDREED (NIN,IREC,A(LMC),NWDS)
*     ----------------------------------
*     MASS RATIO OF EACH MATERIAL/ISOTOPE IN THE CALCULATION DOMAIN:
      AWR(IMT)=A(LMC+MULT)
      NWDS=NWDS+MULT-1
      L3=LMC+NWDS
      L3H=LMCH+NWDS/MULT
      ALLOCATE(TERP(NSUBM*NGRO),TEMP(NSUBM),SIGZ(NSUBM))
      DO 307 ISUBM=1,NSUBM
      TEMP(ISUBM)=A(LMC+MULT+6*(ISUBM-1)+1)
      SIGZ(ISUBM)=A(LMC+MULT+6*(ISUBM-1)+2)
  307 CONTINUE
      CALL LIBTER(NGRO,NSUBM,TEMP,SIGZ,TN(IMT),SN(1,IMT),TERP)
      DEALLOCATE(SIGZ,TEMP)
      L5=0
      IFTOT=0
*----
*  TEMPERATURE AND BACKGROUND LOOP
*----
      DO 600 ISUBM=1,NSUBM
      LOC=LMC+MULT+6*(ISUBM-1)
      TMAT=A(LOC+1)
      SMAT=A(LOC+2)
      LOCS=IA(LOC+6)
      LSUBM1=(ISUBM.EQ.1)
      IF(.NOT.LSUBM1) THEN
         LTERP=.TRUE.
         DO 324 IK=1,NGRO
         LTERP=LTERP.AND.(TERP(NGRO*(ISUBM-1)+IK).EQ.0.0)
  324    CONTINUE
         IF(LTERP) GO TO 600
      ENDIF
*----
*  PROCESS THIS SUBMATERIAL
*----
      LOC=LMC+MULT+6*(ISUBM-1)
      N1DR=IA(LOC+3)
      N1DB=IA(LOC+4)
      N2DB=IA(LOC+5)
      JREC=IREC+LOCS
*----
*  VECTOR CONTROL
*----
      IF(N1DR.EQ.0) GO TO 475
      NWDS=(3+MULT)*N1DR
      IF(L3+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(5).')
      JREC=JREC+1
*     ---------------------------------
      CALL XDREED (NIN,JREC,A(L3),NWDS)
*     ---------------------------------
      NEX1=L3-1+MULT*N1DR
      NEX2=NEX1+N1DR
      NEX3=NEX2+N1DR
      IF(LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,870) HTYPE,HMAT,(HA(L3H+IR-1),IR=1,N1DR)
      ENDIF
*----
*  VECTOR PARTIALS
*----
      IF(LSUBM1) THEN
         IFTOT=0
*        IF NF IS PRESENT, SET IFTOT=1 AND USE NF+NNF+N2NF+N3NF
         DO 325 IR=1,N1DR
         WRITE(HVPS,FORM) HA(L3H-1+IR)
         IF(HVPS.EQ.'NF') IFTOT=1
  325    CONTINUE
         DO 335 KG=1,NGRO
         SFIS(KG)=0.0
         SAVE(KG)=0.0
  335    CONTINUE
      ENDIF
*----
*  LOOP OVER REACTIONS
*----
      IB=0
      DO 470 IR=1,N1DR
      IBLK=IA(NEX1+IR)
      IF(IBLK.GT.IB) THEN
         NWDS=0
*        MANY VECTORS (REACTIONS) ARE STORED IN BLOCK IBLK.
         DO 340 IJ=1,N1DR
         IF(IA(NEX1+IJ).NE.IBLK) GO TO 340
         NWDS=NWDS+IA(NEX3+IJ)-IA(NEX2+IJ)+1
  340    CONTINUE
         ALLOCATE(XS(NWDS))
         JREC=JREC+1
*        ------------------------------
         CALL XDREED (NIN,JREC,XS,NWDS)
*        ------------------------------
         IB=IBLK
         L5=0
      ENDIF
      WRITE(HVPS,FORM) HA(L3H-1+IR)
      NK=IA(NEX3+IR)-IA(NEX2+IR)+1
*----
*  SAVE REQUIRED EXTRA EDIT.
*----
      DO 346 IED=1,NED
      IF(HVPS.EQ.HVECT(IED)) THEN
         IF(LSUBM1) THEN
            DO 341 IK=1,NGRO
            VECT(IK)=0.0
  341       CONTINUE
         ELSE
            CALL LCMGET(KPLIB,HVECT(IED),VECT)
         ENDIF
         DO 345 IK=1,NK
         IF(XS(L5+IK).EQ.0.0) GO TO 345
         JJ=IA(NEX2+IR)+IK-1
         TERPZ=1.0
         IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJ)
         VECT(JJ)=VECT(JJ)+TERPZ*XS(L5+IK)
  345    CONTINUE
         LOGIED(IED,IMT)=.TRUE.
         CALL LCMPUT(KPLIB,HVECT(IED),NGRO,2,VECT)
         GO TO 347
      ENDIF
  346 CONTINUE
*----
*  SAVE MODEL WEIGHT FUNCTIONS
*----
  347 IF((HTYPE.EQ.'NSCAT').AND.(HVPS.EQ.'NWT0').AND.LSUBM1) THEN
         DO 355 IK=1,NK
         JJ=IA(NEX2+IR)+IK-1
         FLUX(JJ,IMT)=XS(L5+IK)
  355    CONTINUE
         GO TO 466
      ENDIF
      IF((HTYPE.EQ.'NTHERM').AND.(HVPS.NE.HINC).AND.
     1   (HVPS.NE.HCOH)) GO TO 466
*----
*  LOOP OVER GROUPS
*----
      DO 440 IK=1,NK
      IF(XS(L5+IK).EQ.0.0) GO TO 440
      JJ=IA(NEX2+IR)+IK-1
      LTHERM=(JJ.GE.NGRO-NTFG(IMT)+1)
      LTIME=(ITIME.EQ.1)
*----
*  INTERPOLATION FACTOR
*----
      TERPZ=1.0
      IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJ)
      IF((SMAT.LT.0.9E10).AND.(ABS(XS(L5+IK)).GT.1.0E-6).AND.
     1 (.NOT.LSUBM1).AND.(HVPS.EQ.'NTOT0')) THEN
         NGF=MIN(NGF,JJ-1)
         NGFR=MAX(NGFR,JJ)
      ENDIF
      IF(ABS(TERPZ).LT.1.0E-3) GO TO 440
      ADD=TERPZ*XS(L5+IK)
*
      IF(HVPS.EQ.'NTOT0') THEN
*        TOTAL XSEC
         TOTAL(JJ,IMT)=TOTAL(JJ,IMT)+ADD
      ELSE IF((.NOT.LSUBM1).AND.(HVPS.EQ.'NFTOT')) THEN
*        FISSION CROSS SECTION
         SIGF(JJ,IMT)=SIGF(JJ,IMT)+ADD*SAVE(JJ)
      ELSE IF(LSUBM1.AND.(HVPS.EQ.'NFTOT')) THEN
         SFIS(JJ)=SFIS(JJ)+ADD
      ELSE IF(LSUBM1.AND.(HVPS.EQ.'NFSLO')) THEN
*        SLOW FISSION
         SIGF(JJ,IMT)=SIGF(JJ,IMT)+ADD
         SAVE(JJ)=SAVE(JJ)+ADD
         IF(IK.EQ.1) SNORM(IMT)=0.0
         SNORM(IMT)=SNORM(IMT)+ADD*FLUX(JJ,IMT)
      ELSE IF(LSUBM1.AND.(HVPS.EQ.'CHIS')) THEN
*        SLOW FISSION
         IF(SNORM(IMT).EQ.0.0) THEN
            WRITE (HSMG,1050) HMAT
            CALL XABORT(HSMG)
         ENDIF
         ADDD=SNORM(IMT)*XS(L5+IK)
         CNORM(IMT)=CNORM(IMT)+ADDD
         CHI(JJ,IMT)=CHI(JJ,IMT)+ADDD
      ELSE IF(LSUBM1.AND.LTIME.AND.(HVPS.EQ.'NUDEL')) THEN
*        DELAYED FISSION
         SIGF(JJ,IMT)=SIGF(JJ,IMT)+ADD*SFIS(JJ)
         SAVE(JJ)=SAVE(JJ)+SFIS(JJ)*ADD
         IF(IK.EQ.1) DNORM(IMT)=0.0
         DNORM(IMT)=DNORM(IMT)+ADD*SFIS(JJ)*FLUX(JJ,IMT)
      ELSE IF(LSUBM1.AND.LTIME.AND.(HVPS.EQ.'CHID')) THEN
*        DELAYED FISSION
         IF(DNORM(IMT).EQ.0.0) THEN
            WRITE (HSMG,1060) HMAT
            CALL XABORT(HSMG)
         ENDIF
         ADDD=DNORM(IMT)*XS(L5+IK)
         CNORM(IMT)=CNORM(IMT)+ADDD
         CHI(JJ,IMT)=CHI(JJ,IMT)+ADDD
      ENDIF
  440 CONTINUE
*
*     END OF REACTION LOOP
  466 L5=L5+NK
      IF(L5.EQ.NWDS) DEALLOCATE(XS)
  470 CONTINUE
*----
*  SCATTERING MATRIX CONTROL
*----
  475 IF(N2DB.EQ.0) GO TO 600
      DO 580 K=1,N2DB
      NWDS=MULT+2+2*NOUTG
      IF(L3+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBTR1: INSUFFICIENT VALUE OF MAXA(6).')
      JREC=JREC+1
*     ---------------------------------
      CALL XDREED (NIN,JREC,A(L3),NWDS)
*     ---------------------------------
      LORD=IA(L3+MULT+1)
      IF(LORD.EQ.0) GO TO 580
      WRITE(HMTX,FORM) HA(L3H)
      LONE=IA(L3+MULT)
      LN=L3+MULT+1
      LG=LN+NOUTG
      IFISN=0
      IF(HTYPE.EQ.'NSCAT'.AND.(HMTX.EQ.'NF'.OR.HMTX.EQ.'NNF'
     1   .OR.HMTX.EQ.'N2NF'.OR.HMTX.EQ.'N3NF')) IFISN=1
      IF(HTYPE.EQ.'NSCAT'.AND.HMTX.EQ.'NFTOT')IFISN=2
*----
*  SCATTERING SUB-BLOCKS
*----
      INC=(NOUTG-1)/NSBLK+1
      DO 570 J=1,NSBLK
      NWDS=0
      DO 480 JJ=(J-1)*INC+1,MIN(J*INC,NOUTG)
      NWDS=NWDS+IA(LN+JJ)
  480 CONTINUE
      IF(NWDS.EQ.0) GO TO 570
      NWDS=NWDS*LORD
      ALLOCATE(XS(NWDS))
      JREC=JREC+1
*     ------------------------------
      CALL XDREED (NIN,JREC,XS,NWDS)
*     ------------------------------
      IF(IFTOT.EQ.1.AND.IFISN.EQ.2) GO TO 560
*----
*  STORE DESIRED CROSS SECTIONS
*----
      IF(HTYPE.EQ.'NTHERM'.AND.HMTX.NE.HINC.AND.
     1   HMTX.NE.HCOH) GO TO 530
      L5=0
*----
*  LOOP OVER SINK, ORDER, SOURCE
*----
      DO 525 JJ=(J-1)*INC+1,MIN(J*INC,NOUTG)
      NP=IA(LN+JJ)
      IF(NP.EQ.0) GO TO 520
      DO 510 IL=1,LORD
      ILNOW=IL+LONE
      IF(ILNOW.GT.NL) GO TO 510
      DO 500 IP=1,NP
      XSNOW=XS(L5+IP+NP*(IL-1))
      IF(XSNOW.EQ.0.) GO TO 500
      JJP=IA(LG+JJ)-IP+1
*----
*  INTERPOLATION FACTOR
*----
      TERPZ=1.0
      IF(.NOT.LSUBM1) TERPZ=TERP(NGRO*(ISUBM-1)+JJP)
      IF(ABS(TERPZ).LT.1.0E-3) GO TO 500
      XSEC=TERPZ*XSNOW
*----
*  CHECK FOR FISSION MATRICES
*----
      IF(IFISN.GT.0) GO TO 490
*----
*  THERMAL CORRECTION TO SCATTERING MATRIX
*----
      IF((HMTX.EQ.'NELAS').AND.(JJP.GE.NGRO-NTFG(IMT)+1)) THEN
         IF(ILNOW.EQ.1) TOTAL(JJP,IMT)=TOTAL(JJP,IMT)-XSEC
         GO TO 500
      ENDIF
      IF(((HMTX.EQ.HINC).OR.(HMTX.EQ.HCOH)).AND.(JJP.LT.
     1 NGRO-NTFG(IMT)+1)) GO TO 500
*----
*  TOTAL SCATTERING MATRIX
*----
*     SCAT(SECONDARY,PRIMARY,ORDER+1)
      SCAT(JJ,JJP,ILNOW)=SCAT(JJ,JJP,ILNOW)+XSEC
*----
*  TOTAL XS AND TOTAL SCATTERING VECTOR
*----
      SIGS(JJP,ILNOW,IMT)=SIGS(JJP,ILNOW,IMT)+XSEC
      IF((ILNOW.EQ.1).AND.(JJP.GE.NGRO-NTFG(IMT)+1)) THEN
         TOTAL(JJP,IMT)=TOTAL(JJP,IMT)+XSEC
      ENDIF
*----
*  FISSION VECTORS
*----
  490 IF(ILNOW.NE.1) GO TO 500
      IF(IFTOT.EQ.1.AND.IFISN.NE.1) GO TO 500
      IF(IFTOT.EQ.0.AND.IFISN.NE.2) GO TO 500
      SIGF(JJP,IMT)=SIGF(JJP,IMT)+XSEC
      CNORM(IMT)=CNORM(IMT)+XSEC*FLUX(JJP,IMT)
      CHI(JJ,IMT)=CHI(JJ,IMT)+XSEC*FLUX(JJP,IMT)
  500 CONTINUE
  510 CONTINUE
  520 L5=L5+NP*LORD
  525 CONTINUE
*----
*  ACCUMULATE FISSION NUBAR
*----
  530 IF(LSUBM1.AND.(HTYPE.EQ.'NSCAT')) THEN
         IF(IFTOT.EQ.1.AND.IFISN.NE.1) GO TO 560
         IF(IFTOT.EQ.0.AND.IFISN.NE.2) GO TO 560
         L5=0
         DO 555 JJ=(J-1)*INC+1,MIN(J*INC,NOUTG)
         NP=IA(LN+JJ)
         IF(NP.EQ.0) GO TO 550
         DO 540 IP=1,NP
         JJP=IA(LG+JJ)-IP+1
         SAVE(JJP)=SAVE(JJP)+XS(L5+IP)
  540    CONTINUE
  550    L5=L5+NP*LORD
  555    CONTINUE
      ENDIF
  560 DEALLOCATE(XS)
  570 CONTINUE
      HGAR(MOD(K-1,18)+1)=HMTX
      IF((K.EQ.1).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,880) HTYPE,HMAT
      ENDIF
      IF((MOD(K-1,18).EQ.17).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,885) (HGAR(I)//' ',I=1,18)
      ELSE IF((K.EQ.N2DB).AND.LSUBM1.AND.(IMPX.GT.4)) THEN
         WRITE (IOUT,885) (HGAR(I)//' ',I=1,MOD(N2DB-1,18)+1)
      ENDIF
  580 CONTINUE
*----
*  SAVE FISSION NU FOR SHIELDING TERMS
*----
      IF(LSUBM1.AND.(HTYPE.EQ.'NSCAT')) THEN
         DO 590 JJ=1,NGRO
         IF(SFIS(JJ).EQ.0) GO TO 590
         SAVE(JJ)=SAVE(JJ)/SFIS(JJ)
  590    CONTINUE
      ENDIF
*----
*  END OF SUBMATERIAL LOOP
*----
  600 CONTINUE
      DEALLOCATE(TERP)
*----
*  SAVE SCATTERING MATRICES ON LCM
*----
      CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGRO,SIGS(1,1,IMT),SCAT,ITYPRO)
*
      GO TO 300
*----
*  END OF MATERIAL AND DATA TYPE LOOPS
*----
  670 CONTINUE
  680 CONTINUE
*----
*  CLOSE MATXS FILE.
*----
      CALL XDRCLS(NIN)
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(37HLIBTR1: UNABLE TO CLOSE LIBRARY FILE ,A8,1H.
     1   )') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
*----
*  CHECK IF ALL NBISO ISOTOPES HAVE BEEN PROCESSED.
*----
      NISOT=0
      DO 700 I=1,NBISO
      IF(MASKI(I)) THEN
         IF((IPR(1,I).EQ.0).OR.(IPR(2,I).EQ.0)) THEN
            WRITE (IOUT,910) (ISONAM(ITC,I),ITC=1,3),NAMFIL
            NISOT=NISOT+1
         ENDIF
      ENDIF
  700 CONTINUE
      IF(NISOT.GT.0) CALL XABORT('LIBTR1: MISSING ISOTOPES')
*----
*  PRINT FINAL FLUX COMPONENTS
*----
      IF(IMPX.GT.6) THEN
         DO 720 IRG=1,NBISO
         IF(MASKI(IRG)) THEN
            SUM=0.0
            DO 710 JJ=1,NGRO
            SUM=SUM+FLUX(JJ,IRG)
  710       CONTINUE
            WRITE(IOUT,927) (ISONAM(ITC,IRG),ITC=1,3),SUM
            WRITE(IOUT,928) (FLUX(I,IRG),I=1,NGRO)
         ENDIF
  720    CONTINUE
      ENDIF
*----
*  PERFORM LIVOLANT-JEANPIERRE NORMALIZATION AND SAVE CROSS SECTION
*  INFORMATION ON LCM.
*----
      DO 830 IM=1,NBISO
      IF(MASKI(IM)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,IM),ITC=1,3)
         KPLIB=IPISO(IM) ! set IM-th isotope
         DO 740 I=1,NGRO
         IF((SN(I,IM).NE.SB(I,IM)).AND.(SN(I,IM).LT.1.0E10)) THEN
            VECT(I)=1.0/(1.0+(TOTAL(I,IM)-SIGS(I,1,IM))*(1.0/SN(I,IM)-
     1      1.0/SB(I,IM)))
         ELSE
            VECT(I)=1.0
         ENDIF
         IF(SN(I,IM).LT.1.0E10) THEN
            FLUX(I,IM)=SN(I,IM)/(SN(I,IM)+TOTAL(I,IM)-SIGS(I,1,IM))/
     1      VECT(I)
         ELSE
            FLUX(I,IM)=1.0
         ENDIF
         TOTAL(I,IM)=TOTAL(I,IM)*VECT(I)
  740    CONTINUE
         IF(IMPX.GT.5) THEN
            WRITE(IOUT,920) HNAMIS
            WRITE(IOUT,928) (VECT(I),I=1,NGRO)
         ENDIF
         CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,TOTAL(1,IM))
         CALL LCMPUT(KPLIB,'NWT0',NGRO,2,FLUX(1,IM))
         CALL XDRLGS(KPLIB,-1,0,0,NL-1,1,NGRO,SIGS(1,1,IM),SCAT,
     1   ITYPRO)
         DO 752 IL=0,NL-1
         DO 751 IG2=1,NGRO
         FACTOR=VECT(IG2)
         SIGS(IG2,IL+1,IM)=SIGS(IG2,IL+1,IM)*FACTOR
         DO 750 IG1=1,NGRO
         SCAT(IG1,IG2,IL+1)=SCAT(IG1,IG2,IL+1)*FACTOR
  750    CONTINUE
  751    CONTINUE
  752    CONTINUE
         CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGRO,SIGS(1,1,IM),SCAT,
     1   ITYPRO)
         DO 780 IED=1,NED
         IF(LOGIED(IED,IM).AND.(HVECT(IED)(:3).NE.'CHI')
     1                    .AND.(HVECT(IED)(:2).NE.'NU')
     2                    .AND.(HVECT(IED).NE.'NTOT0')
     3                    .AND.(HVECT(IED)(:3).NE.'NWT')) THEN
            CALL LCMGET(KPLIB,HVECT(IED),GAR)
            DO 770 I=1,NGRO
            GAR(I)=GAR(I)*VECT(I)
  770       CONTINUE
            CALL LCMPUT(KPLIB,HVECT(IED),NGRO,2,GAR)
         ENDIF
  780    CONTINUE
*
         IF(CNORM(IM).NE.0.0) THEN
*           FISSION SOURCE NORMALIZATION
            DO 790 JJ=1,NGRO
            CHI(JJ,IM)=CHI(JJ,IM)/CNORM(IM)
            SIGF(JJ,IM)=SIGF(JJ,IM)*VECT(JJ)
  790       CONTINUE
            CALL LCMPUT(KPLIB,'NUSIGF',NGRO,2,SIGF(1,IM))
            CALL LCMPUT(KPLIB,'CHI',NGRO,2,CHI(1,IM))
         ENDIF
         CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
         CALL LCMPUT(KPLIB,'AWR',1,2,AWR(IM))
         WRITE(README(:8),'(A8)') HNAMIS(1:8)
         READ(README,'(22A4)') (IHGAR(I),I=1,22)
         CALL LCMPUT(KPLIB,'README',22,3,IHGAR)
      ENDIF
  830 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LOGIED)
      DEALLOCATE(GAR,VECT,FLUX,SCAT,TOTAL,SIGF,SIGS,CHI,SAVE,SFIS,
     1 DNORM,SNORM,CNORM,AWR)
      DEALLOCATE(ITYPRO,IPR)
      RETURN
*
  870 FORMAT(/52H AVAILABLE IDENTIFIERS OF REACTION VECTORS FOR TYPE ,
     1 A6,14H AND MATERIAL ,A6,1H:/(1X,18A7))
  880 FORMAT(/53H AVAILABLE IDENTIFIERS OF REACTION MATRICES FOR TYPE ,
     1 A6,14H AND MATERIAL ,A6,1H:)
  885 FORMAT(1X,18A7)
  890 FORMAT(/32H PROCESSING MATXS LIBRARY NAMED ,A8,1H.)
  910 FORMAT(/27H LIBTR1: MATERIAL/ISOTOPE ',3A4,16H' IS MISSING ON ,
     1 16HMATXS FILE NAME ,A8,1H.)
  920 FORMAT(/40H L-J NORMALIZATION FACTORS FOR MATERIAL ,A12)
  927 FORMAT(/19H FLUX FOR MATERIAL ,3A4,7H   SUM=,1P,E12.5)
  928 FORMAT(1X,1P,10E12.4)
  935 FORMAT(/16H MATXS FILE ID: ,3A6,6H VERS ,I2)
 1050 FORMAT(35HLIBTR1: SNORM MISSING FOR MATERIAL ,A6,1H.)
 1060 FORMAT(35HLIBTR1: DNORM MISSING FOR MATERIAL ,A6,1H.)
      END
