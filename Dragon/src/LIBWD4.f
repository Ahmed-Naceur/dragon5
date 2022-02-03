*DECK LIBWD4
      SUBROUTINE LIBWD4(IPLIB,IPRINT,NAMFIL,NGROUP,NBISO,NL,ISONAM,
     >                  ISONRF,IPISO,ISHINA,TN,SN,SB,MASKI,NGF,NGFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the interpolated microscopic xs read from a
* microscopic xs library in WIMS-D4 format to LCM data structures.
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
* IPLIB   pointer to the internal library.
* IPRINT  print flag.
* NAMFIL  WIMS-D4 library file name.
* NGROUP  number of groups.
* NBISO   number of isotopes.
* NL      number of Legendre scattering order:
*         =1 isotropic;
*         =2 linearly anisotropic;
*         etc.
* ISONAM  local isotope names.
* ISONRF  library isotope names.
* IPISO   pointer array towards microlib isotopes.
* ISHINA  self-shielding isotope names.
* TN      isotope tempterature.
* SN      dilution xs.
* SB      Livolant-Jeanpierre dilution xs.
* MASKI   logical mask for processing isotope.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER    NDPROC
      PARAMETER (NDPROC=10)
      INTEGER    IPRINT,NGROUP,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),
     >           ISHINA(3,NBISO),NGF,NGFR
      CHARACTER  NAMFIL*8,NAMDXS(NDPROC)*6
      LOGICAL    MASKI(NBISO)
      REAL       TN(NBISO),SN(NGROUP,NBISO),SB(NGROUP,NBISO)
*----
* FUNCTIONS
*----
      DOUBLE PRECISION XDRCST
*----
*  INTERNAL PARAMETERS
*----
      INTEGER    IOUT,ITLIB,MAXTEM,MAXDIL,NOTX
      REAL       CONVM
      PARAMETER (IOUT=6,ITLIB=2,MAXTEM=20,MAXDIL=20,NOTX=-1)
      TYPE(C_PTR) KPLIB
      CHARACTER  NAMSBR*6
      PARAMETER (NAMSBR='LIBWD4')
*----
*  LOCAL VARIABLES
*----
      CHARACTER  HNAMIS*12,HSHIR*8
      REAL       TMPT(MAXTEM),DILT(MAXDIL),REST(MAXDIL*MAXTEM),XSCOR(4)
      DOUBLE PRECISION TERP(MAXTEM)
      INTEGER    IP1,NPROC,IUNIT,KDROPN,II,NEL,NGR,NGTHER,MXSCT,IENDF,
     >           ITC,IEL,JEL,JSO,NGX,IG,JC,NRTOT,IELRT,NFIS,NISOR,NSCT,
     >           IT,ILOCX,ILOCY,ILOCS,NRDT,ITXS,IACT,NSRES,IDRES,ILCR,
     >           IXRES,IRES,NTYP,IGF,IGRF,IGR,ITYP,NTMPR,NDILR,ITT,IGRL,
     >           IG1,ISETP1,JG,NBAC,IERR,KDRCLS,IP0
      REAL       ENDR,XX,RIND,XIND,XRS1
*----
*  WIMS-D4 LIBRARY PARAMETERS
*   IUTYPE  type of file = 2 (binary)
*   LRIND   lenght record on da file = 0
*   IACTO   open action = 2 (read only)
*   IACTC   close action = 2 (keep)
*   MAXISO  maximum number of isotopes = 246
*   LPZ     length of Wims parameter array = 8
*   NSETP1  number of p1 scattering sets = 4
*   NPZ     list of main parameters
*   IWISO   id of isotope
*   IDIEL   isotopic id
*   IZ      isotopic charge
*   NF      number fission
*   NR      number resonance
*----
      INTEGER      IUTYPE,LRIND,IACTO,IACTC,MAXISO,LPZ,NSETP1
      PARAMETER   (IUTYPE=2,LRIND=0,IACTO=2,IACTC=1,MAXISO=246,
     >             LPZ=8,NSETP1=4)
      CHARACTER    CWISO(MAXISO)*8,FMT*6
      INTEGER      NPZ(LPZ),IWISO(MAXISO),IDIEL,IZ,NFIEL,
     >             NF(MAXISO),NTMP,NRIEL,NR(MAXISO),IDTEMP(2),
     >             NBATOM(NSETP1)
      REAL         AWR
      INTEGER      IPRLOC
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO,ISORD,NTM,NDI
      REAL, ALLOCATABLE, DIMENSION(:) :: DELTA,XSSCMP,AW,ENER,TMPXS,
     > TMPSC,RID,RTMP,RDIL,RESI,RRI,RIT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSREC,XSOUT,GAR,DSIGPL
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  DATA
*----
      SAVE         NBATOM,NAMDXS
      DATA         NBATOM
     >            /1,2,16,12/
      DATA         NAMDXS
     >            /'NTOT0 ','TRANC ','NUSIGF','NFTOT ','CHI   ',
     >             'NU    ','NG    ','N2N   ','NGOLD ','NWT0  '/
*----
*  SCRATCH STORAGE ALLOCATION
*   ITYPRO  cross section processed
*   DELTA   lethargy
*   XSREC   general xs vector
*   SCAT    complete scattering matrix SCAT(JG,IG) (from IG to JG)
*   XSSCMP  compress scattering for transfer
*   XSOUT   self shielding parameter
*   ISORD   local isotope flag
*   AW      isotope atomic weight
*   GAR     intermediate xs vector:
*           GAR(I,1): library fission spectrum;
*           GAR(I,2): potential scattering xs;
*           GAR(I,3): transport xs;
*           GAR(I,4): absorption xs
*----
      ALLOCATE(ITYPRO(NL),ISORD(NBISO))
      ALLOCATE(DELTA(NGROUP),XSREC(NGROUP,NDPROC+NL),
     >         SCAT(NGROUP,NGROUP,NL),XSSCMP(NGROUP*(NGROUP+2)),
     >         XSOUT(NGROUP,7),AW(NBISO),GAR(NGROUP,4))
*----
*  OPEN WIMS-D4 LIBRARY
*  READ GENERAL DIMENSIONING
*----
      IPRLOC=0
      IF(ABS(IPRINT) .GE. 100) IPRLOC=100
      CONVM=REAL(XDRCST('Neutron mass','amu'))
      IP0=NDPROC+1
      IP1=NDPROC+2
      NPROC=NDPROC+NL
      IUNIT=KDROPN(NAMFIL,IACTO,IUTYPE,LRIND)
      IF(IUNIT.LE.0) CALL XABORT(NAMSBR//': WIMS-D4 LIBRARY '//
     >    NAMFIL//' CANNOT BE OPENED FOR MIXS')
      IF(ABS(IPRINT) .GE. 5) THEN
        WRITE(IOUT,6000) NAMSBR,NAMFIL
      ENDIF
      READ(IUNIT) (NPZ(II),II=1,LPZ)
      IF(NPZ(2).NE.NGROUP) THEN
        WRITE(IOUT,9001) NGROUP,NPZ(2)
        CALL XABORT(NAMSBR//': INVALID NUMBER OF GROUPS')
      ENDIF
      NEL=NPZ(1)
      NGF=NPZ(4)
      NGR=NPZ(5)
      NGTHER=NPZ(6)
      NGFR=NGF+NGR
      MXSCT=NGROUP*(NGROUP+2)
      IF(NGFR+NGTHER.NE.NGROUP) THEN
        WRITE(IOUT,9001) NGROUP,NGFR+NGTHER
        CALL XABORT(NAMSBR//': INVALID NUMBER OF GROUPS')
      ENDIF
      IF(NEL.GT.MAXISO) THEN
        WRITE(IOUT,9002) MAXISO,NEL
        CALL XABORT(NAMSBR//': INVALID NUMBER OF ISOTOPES')
      ENDIF
      IENDF=0
      ALLOCATE(DSIGPL(NGR,NEL))
*----
*  READ ISOTOPE ID NUMBER AND CREATE EQUIVALENT ISOTOPE NAME
*  SCAN TO ASSOCIATE WIMS ISOTOPE NUMBER WITH DRAGON ISOTOPE NUMBER
*  VERIFY IF ALL ISOTOPES REQUIRED ARE PRESENT
*----
      READ(IUNIT) (IWISO(ITC),ITC=1,NEL)
      CALL XDISET(ISORD,NBISO,0)
      DO 100 IEL=1,NEL
        CWISO(IEL)='        '
        IF     (IWISO(IEL).LT.10) THEN
          WRITE(CWISO(IEL),'(I1)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100) THEN
          WRITE(CWISO(IEL),'(I2)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.1000) THEN
          WRITE(CWISO(IEL),'(I3)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.10000) THEN
          WRITE(CWISO(IEL),'(I4)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100000) THEN
          WRITE(CWISO(IEL),'(I5)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.1000000) THEN
          WRITE(CWISO(IEL),'(I6)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.10000000) THEN
          WRITE(CWISO(IEL),'(I7)') IWISO(IEL)
        ELSE IF(IWISO(IEL).LT.100000000) THEN
          WRITE(CWISO(IEL),'(I8)') IWISO(IEL)
        ENDIF
        READ(CWISO(IEL),'(2A4)') (IDTEMP(ITC),ITC=1,2)
        DO 101 JSO=1,NBISO
          IF(MASKI(JSO)) THEN
            IF(ISONRF(1,JSO).EQ.IDTEMP(1).AND.
     >         ISONRF(2,JSO).EQ.IDTEMP(2)) ISORD(JSO)=IEL
          ENDIF
 101    CONTINUE
 100  CONTINUE
      DO 102 JSO=1,NBISO
        IF(MASKI(JSO).AND.(ISORD(JSO).EQ.0)) THEN
          WRITE(IOUT,9003) (ISONRF(ITC,JSO),ITC=1,3),NAMFIL
          CALL XABORT(NAMSBR//': MISSING ISOTOPE')
        ENDIF
 102  CONTINUE
*----
*  READ GROUP STRUCTURE
*----
      ALLOCATE(ENER(NGROUP+1))
      READ(IUNIT) (ENER(ITC),ITC=1,NGROUP+1)
      IF(ENER(NGROUP+1).EQ.0.0) ENER(NGROUP+1)=1.0E-5
      CALL LCMPUT(IPLIB,'ENERGY',NGROUP+1,2,ENER)
      NGX=0
      DO 103 IG=1,NGROUP
        IF(NGX.EQ.0.AND.ENER(IG+1).LT.4.0) NGX=IG-1
        DELTA(IG)=LOG(ENER(IG)/ENER(IG+1))
 103  CONTINUE
      DEALLOCATE(ENER)
      CALL LCMPUT(IPLIB,'DELTAU',NGROUP,2,DELTA)
*----
*  RECOVER FISSION SPECTRUM
*  INITIALIZE FAST AND THERMAL GROUPS XPO AND GOLD
*  ONLY RESONANCE GROUPS FOR THESE CROSS SECTIONS
*  INITIALIZE XN2N
*  NO SUCH REACTION FOR THIS LIBRARY
*----
      CALL XDRSET(GAR(1,1),NGROUP,0.0)
      READ(IUNIT) (GAR(ITC,1),ITC=1,NPZ(3))
*----
*  READ DEPLETION CHAIN
*----
      DO 120 IEL=1,NEL
        READ(IUNIT) JC
 120  CONTINUE
      READ(IUNIT) ENDR
*----
*  ALLOCATE MEMORY FOR TEMPERATURE DEPENDENT XS
*  AND FOR RESONANCE CALCULATION
*----
      ALLOCATE(TMPXS(NGROUP*5*MAXTEM),TMPSC(NGROUP*NGROUP*MAXTEM))
*----
* READ FILE
* CROSS SECTION ARE SAVED ONLY IF ISOTOPE IS USED
*----
      CALL XDRSET(AW,NBISO,0.0)
      NRTOT=0
      DO 130 IELRT=1,NEL
        READ(IUNIT) IDIEL,AWR,IZ,NFIEL,NTMP,NRIEL
        IF(NRIEL.GT.0) THEN
          NRTOT=NRTOT+NRIEL
        ENDIF
        IF(NTMP.GT.MAXTEM) THEN
          CALL XABORT(NAMSBR//': INVALID MAXTEM FOR P0.')
        ENDIF
*----
*  LOCATE ISOTOPE IN LIST OF LIBRARY ISOTOPES IN THE CASE
*  WHERE LIBRARY IS NOT COMPLETE OR THE ORDER OF ISOTOPE
*  STORED IS DIFFERENT FROM THAT OF THE ISOTOPE NAMES
*----
        IEL=0
        DO 140 JEL=1,NEL
          IF(IDIEL.EQ.IWISO(JEL)) THEN
            IEL=JEL
            NF(IEL)=NFIEL
            NFIS=0
            IF(NF(IEL).GT.1) NFIS=1
            NR(IEL)=NRIEL
            GO TO 145
          ENDIF
 140    CONTINUE
        CALL XABORT(NAMSBR//': WIMSD4 LIBRARY INCOMPLETE')
 145    CONTINUE
        NISOR=0
*----
*  SCAN TO SEE IF ISOTOPE IS REQUIRED
*----
        DO 150 JSO=1,NBISO
          IF(MASKI(JSO).AND.(ISORD(JSO).EQ.IEL)) THEN
             NISOR=1
             GO TO 155
          ENDIF
 150    CONTINUE
 155    CONTINUE
        IF(NISOR.EQ.0) THEN
*----
*  ISOTOPE NOT REQUIRED/SKIP RECORDS
*----
          READ(IUNIT) XX
          IF(NF(IEL).GT.1) READ(IUNIT) XX
          READ(IUNIT) NSCT
          IF(NTMP.GT.0) THEN
            READ(IUNIT) XX
            DO 160 IT=1,NTMP
              READ(IUNIT) XX
              IF(NF(IEL).GT.1) THEN
                READ(IUNIT) XX
              ENDIF
              READ(IUNIT) NSCT
 160        CONTINUE
          ENDIF
        ELSE
*----
*  ISOTOPE REQUIRED READ FAST AND/OR RESONANCE XS
*----
          CALL XDRSET(XSREC(1,1),NGROUP*NPROC,0.0)
          CALL XDRSET(XSREC(1,9),NGROUP,1.0)
          READ(IUNIT) (GAR(NGF+II,2),II=1,NGR),
     >                (XX,II=1,NGR),
     >                (GAR(II,3),II=1,NGFR),
     >                (GAR(II,4),II=1,NGFR),
     >                (XX,II=1,NGR),
     >                (XSREC(NGF+II,9),II=1,NGR)
          CALL XDRSET(DSIGPL(1,IEL),NGR,0.0)
          DO 180 IG=NGF+1,NGFR
            DSIGPL(IG-NGF,IEL)=GAR(IG,2)*XSREC(IG,9)
 180      CONTINUE
          IF(NF(IEL).GT.1) THEN
            READ(IUNIT) (XSREC(II,3),II=1,NGFR),
     >                  (XSREC(II,4),II=1,NGFR)
            DO 185 IG=1,NGROUP
              XSREC(IG,5)=GAR(IG,1)
 185        CONTINUE
          ENDIF
*----
*  READ AND DECOMPRESS P0 SCATTERING CROSS SECTIONS
*  COMPUTE P0 SCATTERING OUT OF GROUP
*----
          READ(IUNIT) NSCT,(XSSCMP(II),II=1,NSCT)
          CALL LIBWSC(NGROUP,1,NGFR,NSCT,XSSCMP,SCAT,XSREC(1,IP0))
*----
*  THERMAL XS
*----
          IF(NTMP.EQ.1) THEN
            READ(IUNIT) XX
            READ(IUNIT) (GAR(NGFR+II,3),II=1,NGTHER),
     >                  (GAR(NGFR+II,4),II=1,NGTHER)
            IF(NF(IEL).GT.1) THEN
              READ(IUNIT) (XSREC(NGFR+II,3),II=1,NGTHER),
     >                    (XSREC(NGFR+II,4),II=1,NGTHER)
            ENDIF
            READ(IUNIT) NSCT,(XSSCMP(II),II=1,NSCT)
*----
*  READ AND DECOMPRESS P0 SCATTERING CROSS SECTIONS
*  COMPUTE P0 SCATTERING OUT OF GROUP
*----
            CALL LIBWSC(NGROUP,NGFR+1,NGROUP,NSCT,XSSCMP,
     >                  SCAT,XSREC(1,IP0))
*----
*  SAVE INFORMATION FOR ISOTOPES WITHOUT SELF SHIELDING DATA
*----
            DO 200 JSO=1,NBISO
              IF(MASKI(JSO).AND.(ISORD(JSO).EQ.IEL)) THEN
                WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,JSO),ITC=1,3)
                IF(ABS(IPRINT) .GE. 5) THEN
                  WRITE(IOUT,6001) HNAMIS
                  IF(ABS(IPRINT) .GE. 100) THEN
                    WRITE(IOUT,6200) TN(JSO)
                  ENDIF
                ENDIF
                AW(JSO)=AWR/CONVM
*----
*  BUILT TOTAL CROSS SECTION FROM INFORMATION IN XSNG WHICH IS
*  CURRENTLY ABSORPTION AND SIGS WHICH IS TOTAL SCATTERING
*  OUT OF GROUP
*  COMPUTE REAL NG CROSS SECTION WHICH IS
*  CURRENT NG (ABSORPTION)-FISSION-N2N
*----
                DO 201 IG=1,NGROUP
                  XSREC(IG,1)=GAR(IG,4)+XSREC(IG,IP0)
                  IF(NF(IEL).GT.1) THEN
                    XSREC(IG,7)=GAR(IG,4)-XSREC(IG,8)-XSREC(IG,4)
                  ELSE
                    XSREC(IG,7)=GAR(IG,4)-XSREC(IG,8)
                  ENDIF
                  IF(XSREC(IG,4).NE.0) THEN
                    XSREC(IG,6)=XSREC(IG,3)/XSREC(IG,4)
                  ELSE
                    XSREC(IG,6)=0.0
                  ENDIF
 201            CONTINUE
                IF(IENDF.LT.2) THEN
*----
*  COMPUTE TRANSPORT CORRECTION AND STORE IN TRAN
*----
                  DO 202 IG=1,NGROUP
                    XSREC(IG,2)=XSREC(IG,1)-GAR(IG,3)
 202              CONTINUE
                ELSE
                  DO 203 IG=1,NGROUP
                    XSREC(IG,2)=GAR(IG,3)
 203              CONTINUE
                ENDIF
*----
*  SAVE ISOTOPE INFORMATION
*----
                KPLIB=IPISO(JSO) ! set JSO-th isotope
                CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
                CALL LCMPUT(KPLIB,'AWR',1,2,AW(JSO))
                CALL XDRLGS(KPLIB,1,IPRLOC,0,0,1,NGROUP,
     >                      XSREC(1,IP0),SCAT,ITYPRO)
                CALL XDRLXS(KPLIB,1,IPRLOC,NDPROC,NAMDXS,1,NGROUP,XSREC)
              ENDIF
 200        CONTINUE
          ELSE IF(NTMP.GT.1) THEN
*----
*  READ TEMPERATURE DEPENDENT XS
*----
            READ(IUNIT) (TMPT(II),II=1,NTMP)
            ILOCX=0
            ILOCY=NGFR
            ILOCS=0
            NRDT=NGTHER-1
            DO 210 IT=1,NTMP
              READ(IUNIT) (TMPXS(ILOCY+II+1),II=0,NRDT),
     >                    (TMPXS(ILOCY+II+NGROUP+1),II=0,NRDT)
              IF(NF(IEL).GT.1) THEN
                READ(IUNIT) (TMPXS(ILOCY+II+2*NGROUP+1),II=0,NRDT),
     >                      (TMPXS(ILOCY+II+3*NGROUP+1),II=0,NRDT)
              ENDIF
              READ(IUNIT) NSCT,(XSSCMP(II),II=1,NSCT)
*----
*  READ AND DECOMPRESS P0 SCATTERING CROSS SECTIONS
*  COMPUTE P0 SCATTERING OUT OF GROUP
*  COMPUTE TOTAL XS
*----
              CALL LIBWSC(NGROUP,NGFR+1,NGROUP,NSCT,XSSCMP,
     >                    TMPSC(ILOCS+1),TMPXS(ILOCX+4*NGROUP+1))
              ILOCX=ILOCX+5*NGROUP
              ILOCY=ILOCY+5*NGROUP
              ILOCS=ILOCS+NGROUP*NGROUP
 210        CONTINUE
*----
*  SAVE INFORMATION FOR ISOTOPES
*  NO SELF SHIELDING
*----
            DO 220 JSO=1,NBISO
              IF(MASKI(JSO).AND.(ISORD(JSO).EQ.IEL)) THEN
                WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,JSO),ITC=1,3)
                IF(ABS(IPRINT) .GE. 5)  WRITE(IOUT,6001) HNAMIS
                AW(JSO)=AWR/CONVM
*----
*  FIND TEMPERATURE INTERPOLATION COEFFICIENTS
*  INTERPOLATE IN TEMPERATURE
*----
                CALL LIBLEX(NTMP,TN(JSO),TMPT,NOTX,TERP)
                IF(ABS(IPRINT) .GE. 100) THEN
                  WRITE(IOUT,6201) TN(JSO)
                  WRITE(IOUT,6202) (TMPT(ITC),ITC=1,NTMP)
                  WRITE(IOUT,6203) (TERP(ITC),ITC=1,NTMP)
                ENDIF
                ITXS=1
                IACT=1
                CALL LIBWTE(IACT,ITXS,NGROUP,NGTHER,NTMP,NF(IEL),TERP,
     >                      SCAT,XSREC(1,IP0),GAR(1,4),XSREC(1,3),
     >                      XSREC(1,4),GAR(1,3),TMPXS,TMPSC)
*----
*  BUILT TOTAL CROSS SECTION FROM INFORMATION IN XSNG WHICH IS
*  CURRENTLY ABSORPTION AND SIGS WHICH IS TOTAL SCATTERING
*  OUT OF GROUP
*  COMPUTE REAL NG CROSS SECTION WHICH IS
*  CURRENT NG (ABSORPTION)-FISSION-N2N
*----
                DO 221 IG=1,NGROUP
                  XSREC(IG,1)=GAR(IG,4)+XSREC(IG,IP0)
                  IF(NF(IEL).GT.1) THEN
                    XSREC(IG,7)=GAR(IG,4)-XSREC(IG,8)-XSREC(IG,4)
                  ELSE
                    XSREC(IG,7)=GAR(IG,4)-XSREC(IG,8)
                  ENDIF
                  IF(XSREC(IG,4).NE.0) THEN
                    XSREC(IG,6)=XSREC(IG,3)/XSREC(IG,4)
                  ELSE
                    XSREC(IG,6)=0.0
                  ENDIF
 221            CONTINUE
                IF(IENDF.LT.2) THEN
*----
*  COMPUTE TRANSPORT CORRECTION AND STORE IN TRAN
*----
                  DO 222 IG=1,NGROUP
                    XSREC(IG,2)=XSREC(IG,1)-GAR(IG,3)
 222              CONTINUE
                ELSE
                  DO 223 IG=1,NGROUP
                    XSREC(IG,2)=GAR(IG,3)
 223              CONTINUE
                ENDIF
*----
*  SAVE ISOTOPE INFORMATION
*----
                KPLIB=IPISO(JSO) ! set JSO-th isotope
                CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
                CALL LCMPUT(KPLIB,'AWR',1,2,AW(JSO))
                CALL XDRLGS(KPLIB,1,IPRLOC,0,0,1,NGROUP,
     >                      XSREC(1,IP0),SCAT,ITYPRO)
                CALL XDRLXS(KPLIB,1,IPRLOC,NDPROC,NAMDXS,1,NGROUP,XSREC)
              ENDIF
 220        CONTINUE
          ENDIF
        ENDIF
        READ(IUNIT) ENDR
 130  CONTINUE
*----
*  RELEASE MEMORY FOR TEMPERATURE DEPENDENT XS
*----
      DEALLOCATE(TMPSC,TMPXS)
*----
*  ALLOCATE MEMORY FOR RESONANCE READ
*  READ ALL GROUP AND ALL RESONANCES
*----
      NTYP=2
      ALLOCATE(NTM(NTYP*NRTOT*NGR),NDI(NTYP*NRTOT*NGR))
      ALLOCATE(RID(NRTOT),RTMP(MAXTEM*NTYP*NRTOT*NGR),
     > RDIL(MAXDIL*NTYP*NRTOT*NGR),RESI(MAXDIL*MAXTEM*NTYP*NRTOT*NGR))
      CALL XDRSET(RID,NRTOT,0.0)
      CALL XDISET(NTM,NTYP*NRTOT*NGR,0)
      CALL XDISET(NDI,NTYP*NRTOT*NGR,0)
      CALL XDRSET(RTMP,MAXTEM*NTYP*NRTOT*NGR,0.0)
      CALL XDRSET(RDIL,MAXDIL*NTYP*NRTOT*NGR,0.0)
      CALL XDRSET(RESI,MAXDIL*MAXTEM*NTYP*NRTOT*NGR,0.0)
      CALL LIBWRG(IUNIT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,NSRES,RID,NTM,
     >            NDI,RTMP,RDIL,RESI)
*----
*  ALLOCATE MEMORY FOR RESONANCE PROCESSING
*----
      ALLOCATE(RRI(MAXDIL*MAXTEM*2),RIT(MAXDIL))
*----
*   PROCESS RESONANCES
*----
      IF(ABS(IPRINT) .GE. 5) WRITE(IOUT,6010)
      DO 230 JSO=1,NBISO
        IF(.NOT.MASKI(JSO)) GO TO 235
        IEL=ISORD(JSO)
        IF(IEL.EQ.0) CALL XABORT(NAMSBR//': INVALID VALUE OF ISORD')
        IF(NR(IEL).EQ.0) GO TO 235
        NFIS=0
        IF(NF(IEL).GT.1) NFIS=1
        WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,JSO),ITC=1,3)
        KPLIB=IPISO(JSO) ! set JSO-th isotope
        WRITE(HSHIR,'(2A4)') (ISHINA(ITC,JSO),ITC=1,2)
        IDRES=INDEX(HSHIR,'.')
        IF(IDRES.GT.0) THEN
          WRITE(FMT,'(2H(F,I1,3H.1))') IDRES+1
          READ(HSHIR,FMT) RIND
        ELSE
          RIND=FLOAT(IWISO(IEL))
        ENDIF
*----
*  IDENTIFY RESONANCE SET
*  DEFAULF IS RESONNANCE ID SPECIFIED OR FIRST SET ENCOUNTERED
*----
        ILCR=0
        DO 231 IXRES=1,NSRES
          XIND=RID(ILCR+1)
          IF(IDRES.EQ.0) THEN
            XRS1=FLOAT(INT((XIND+0.01)*10.)-INT(XIND+0.01)*10)/10.
            XRS1=ABS(XIND-XRS1-RIND)
          ELSE
            XRS1=ABS(XIND-RIND)
          ENDIF
          IF(XRS1.LE.0.01) THEN
            IRES=IXRES
            GO TO 236
          ENDIF
          ILCR=ILCR+1
 231    CONTINUE
*----
*  START MODIFICATION: 98/05/05 (G.M.)
*  1) NO SPECIFIC RESONNANCE ID SPECIFIED AND
*     NO RESONNANCE ID FOUND WITH ADEQUATE IWISO
*     EVEN IF NR(IEL) > 0
*     ASSUME NO RESONANCE PRESENT
*     REQUIRED FOR ANL 69 GROUPS WIMSD4 LIBRARY
*     ISOTOPE '10' HAS NRES=1 BUT ID DIFFERENT
*     FROM 10.X WHERE x IS AN INTEGER
*  2) CORRECT ERROR IN LOOP INDEX FOR
*     WRITE(IOUT,9004)
*----
        IF(IDRES.EQ.0) GO TO 235
        WRITE(IOUT,9004) (ISONAM(ITC,JSO),ITC=1,3),RIND
*----
*  END MODIFICATION: G.M. (98/05/05)
*----
        CALL XABORT(NAMSBR//': UNABLE TO IDENTIFY RESONANCE SET '//
     >              'FOR THIS ISOTOPE')
 236    CONTINUE
*----
*  THIS ISOTOPE NEEDS TO BE CORRECTED FOR SELF SHIELDING
*  FIRST READ UNCORRECTED CROSS SECTIONS
*----
        NTYP=2
        XSCOR(1)=0.0
        XSCOR(2)=0.0
        XSCOR(3)=0.0
        XSCOR(4)=0.0
        IF(ABS(IPRINT) .GE. 5) WRITE(IOUT,6011) HNAMIS,XIND,TN(JSO)
        CALL XDRLGS(KPLIB,-1,0,0,0,1,NGROUP,XSREC(1,IP0),SCAT, ITYPRO)
        CALL XDRLXS(KPLIB,-1,0,NDPROC,NAMDXS,1,NGROUP,XSREC)
*----
*  SCAN RESONAMCE GROUPS AND CORRECT CROSS SECTIONS
*----
        DO 232 IGF=1,NGFR
          XSOUT(IGF,3)=XSREC(IGF,IP0)
          XSOUT(IGF,4)=1.0
          XSOUT(IGF,5)=1.0
 232    CONTINUE
        IGRF=NGF
        DO 240 IGR=1,NGR
          IGRF=IGRF+1
*----
*  PREPARE VECTORS FOR SELF SHIELDING
*----
          IF(ABS(IPRINT) .GE. 100) THEN
            WRITE(IOUT,6004) IGRF,SN(IGRF,JSO),DSIGPL(IGR,IEL)
          ENDIF
          DO 250 ITYP=1,NTYP
            IF((ITYP.EQ.2).AND.(NF(IEL).NE.3)) GO TO 250
            CALL LIBWRP(IPRINT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,IGR,IRES,
     >                  ITYP,DSIGPL(IGR,IEL),NTM,NDI,RTMP,RDIL,RESI,
     >                  NTMPR,NDILR,TMPT,DILT,REST)
            IF(NDILR.GT.0.AND.NTMPR.GT.0) THEN
              CALL LIBWRI(NTMPR,NDILR,TN(JSO),SN(IGRF,JSO),TMPT,DILT,
     >                    REST,RIT,XSOUT(IGRF,ITYP),XSCOR(ITYP))
              IF(ABS(IPRINT) .GE. 100) THEN
                IF(ITYP.EQ.1) THEN
                  WRITE(IOUT,6002) 'absorption  '
                ELSE IF(ITYP.EQ.2) THEN
                  WRITE(IOUT,6002) 'fission     '
                ELSE IF(ITYP.EQ.3) THEN
                  WRITE(IOUT,6002) 'scattering  '
                ENDIF
                WRITE(IOUT,6003) (REST(ITT),ITT=1,NTMPR*NDILR)
                IF(ITYP.EQ.1) THEN
                  WRITE(IOUT,6005) XSOUT(IGRF,ITYP)
                ELSE IF(ITYP.EQ.2) THEN
                  WRITE(IOUT,6006) XSOUT(IGRF,ITYP)
                ELSE IF(ITYP.EQ.3) THEN
                  WRITE(IOUT,6007) XSOUT(IGRF,ITYP)
                ENDIF
              ENDIF
            ENDIF
 250      CONTINUE
 240    CONTINUE
*----
*  CORRECT CROSS SECTIONS FOR ALL RESONANCE GROUPS
*----
        IGRF=NGF+1
        IGRL=NGF+NGR
        IF(NF(IEL).NE.3) NTYP=1
        CALL LIBWRE(NTYP,IPRINT,ITLIB,NGROUP,1,IGRF,IGRL,NGR,
     >              SCAT,XSREC(1,IP0),XSREC(1,1),XSREC(1,7),
     >              XSREC(1,3),XSREC(1,4),XSREC(1,6),
     >              DELTA,SN(1,JSO),SB(1,JSO),XSOUT,XSCOR,
     >              DSIGPL(1,IEL))
*----
*  PRINT CROSS SECTIONS IF REQUIRED
*----
        IF(ABS(IPRINT) .GE. 5) THEN
          WRITE(IOUT,6100) HNAMIS
          DO 233 IG1=NGF+1,NGFR
            WRITE(IOUT,6101) IG1,SN(IG1,JSO),SB(IG1,JSO),
     >                       XSOUT(IG1,4),XSREC(IG1,1),
     >                       XSREC(IG1,IP0),XSREC(IG1,3),XSREC(IG1,9)
 233      CONTINUE
        ENDIF
*----
*  SET NWT0 THE RESONANCE FLUX WEIGHTING
*----
        CALL XDRSET(XSREC(1,10),NGROUP,1.0)
        DO 234 IG1=NGF+1,NGFR
          XSREC(IG1,10)=XSOUT(IG1,4)
 234    CONTINUE
*----
*  SAVE SELF-SHIELDED XS
*----
        CALL XDRLGS(KPLIB,1,0,0,0,1,NGROUP,XSREC(1,IP0),SCAT, ITYPRO)
        CALL XDRLXS(KPLIB,1,0,NDPROC,NAMDXS,1,NGROUP,XSREC)
 235    CONTINUE
 230  CONTINUE
*----
*  RELEASE MEMORY FOR RESONANCE PROCESSING
*----
      DEALLOCATE(RIT,RRI,RID)
*----
*  RELEASE MEMORY FOR RESONANCE READ
*----
      DEALLOCATE(RESI,RDIL,RTMP)
      DEALLOCATE(NDI,NTM)
*----
*  P1 SCATTERING
*----
      IF(NL.GE.2) THEN
        DO 260 ISETP1=1,NSETP1
          DO 270 IG=1,NGROUP
            XSREC(IG,IP1)=0.0
            READ(IUNIT,END=300) (SCAT(JG,IG,2),JG=1,NGROUP)
            DO 280 JG=1,NGROUP
              XSREC(IG,IP1)=XSREC(IG,IP1)+SCAT(JG,IG,2)
 280        CONTINUE
 270      CONTINUE
          DO 290 JSO=1,NBISO
            NBAC=INT(AW(JSO)*CONVM+0.01)
            IF(NBAC.EQ.NBATOM(ISETP1)) THEN
*----
*  CLASSIFY BY ATOMIC NUMBER
*----
              WRITE(HNAMIS,'(3A4)') (ISONAM(ITC,JSO),ITC=1,3)
              IF(IPRINT.GT.5) WRITE(IOUT,6300) HNAMIS
              KPLIB=IPISO(JSO) ! set JSO-th isotope
              CALL XDRLGS(KPLIB,-1,0,0,0,1,NGROUP,XSREC(1,IP0),
     >                    SCAT,ITYPRO)
              CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGROUP,XSREC(1,IP0),
     >                    SCAT,ITYPRO)
            ENDIF
 290      CONTINUE
 260    CONTINUE
      ENDIF
 300  IERR=KDRCLS(IUNIT,IACTC)
      IF(IERR.LT.0)
     >  CALL XABORT(NAMSBR//': WIMS-D4 LIBRARY '//
     >  NAMFIL//' CANNOT BE CLOSED')
      IF(ABS(IPRINT) .GE. 5) THEN
        WRITE(IOUT,6009) NAMSBR
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DSIGPL)
      DEALLOCATE(GAR,AW,XSOUT,XSSCMP,SCAT,XSREC,DELTA)
      DEALLOCATE(ISORD,ITYPRO)
*----
*  RETURN
*----
      RETURN
*----
*  FORMAT
*----
 9001 FORMAT(/' NUMBER OF GROUPS SPECIFIED :',I10/
     >        ' NUMBER OF GROUPS IN LIBRARY :',I10)
 9002 FORMAT(/' MAXIMUM NUMBER OF ISOTOPE SPECIFIED :',I10/
     >        '        NUMBER OF ISOTOPE IN LIBRARY :',I10)
 9003 FORMAT(/' LIBWD4: MATERIAL/ISOTOPE ',3A4,
     >        ' IS MISSING ON WIMS-D4 FILE ',A8)
 9004 FORMAT(/' LIBWD4: FOR ISOTOPE ',3A4,
     >        ' SELF-SHIELDING ISOTOPE ',F8.1,' NOT AVAILABLE')
 6000 FORMAT('(* Output from --',A6,'-- follows '//
     >       ' READING WIMS-D4 LIBRARY NAME ',A8)
 6001 FORMAT('   PROCESSING ISOTOPE/MATERIAL = ',A12)
 6002 FORMAT('   Resonance integral tabulation for ',A12)
 6003 FORMAT(1P,5E15.7)
 6004 FORMAT('   Processing GROUP = ', I10,'   at  dilutions = ',
     >       1P,2E15.7)
 6005 FORMAT('   Interpolated absorption rate   = ',1P,E15.7)
 6006 FORMAT('   Interpolated fission rate      = ',1P,E15.7)
 6007 FORMAT('   Interpolated scattering rate   = ',1P,E15.7)
 6009 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('   RESONANCE IDENTIFICATION')
 6011 FORMAT('   ISOTOPE ID = ',A12,' RESONANCE ID = ',F8.1,
     >       ' at temperature  = ',F10.5)
 6100 FORMAT('   SELF SHIELDING PROPERTIES FOR ISOTOPE =',A12/
     > 5X,'GROUP',10X,'DILUT',13X,'SB',11X,'NPHI',10X,'NTOT0',
     > 11X,'SIGS',9X,'NUSIGF',10X,'NGOLD')
 6101 FORMAT(5X,I5,1P,8E15.5)
 6200 FORMAT(' TEMPERATURE = ',F10.5,10X,
     >  ' CROSS SECTION TABULATED AT A SINGLE TEMPERATURE')
 6201 FORMAT(' TEMPERATURE = ',F10.5,10X,
     >  ' CROSS SECTION TABULATED AT MULTIPLE TEMPERATURES')
 6202 FORMAT(' TABULATION TEMPERATURES= ',/(5F15.5))
 6203 FORMAT(' INTERPOLATION FACTORS  = ',1P,/(5E15.5))
 6300 FORMAT(' P1 SCATTERING CROSS SECTION FOUND FOR =',A12)
      END
